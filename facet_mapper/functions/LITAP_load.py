#imports
from dbfread import DBF
import pandas as pd
import numpy as np
import copy
import string
from datar import dplyr,tidyr
from datar.all import f
from itertools import compress
import math
import datar

#import from other files
# from LITAP_functions import *   
# from LITAP_utils import *

import LITAP_functions   
import LITAP_utils


#LITAP_load.R functions
def load_dem(file,type_="elev"):
    db = pd.DataFrame(iter(DBF(file, load=True)))
    db.columns = map(str.lower, db.columns)
    #TODO check name function in R
    return db

def db_format(db,nrow,ncol,missing_value=-9999,verbose=True):
    
    if (nrow*ncol!=len(db)):
        print("TON POULO")
    db["seqno"] = range(1, len(db) + 1)
    
    #division may cause problems if numbers of nrow, ncol and len(db)  are not optimal
    db["row"] = sorted(list(range(1,nrow +1)) * (len(db)//nrow))
    db["col"] = list(range(1,ncol +1)) * (len(db)//ncol)
#     print(db)
    db["missing"] = db["elev"]==missing_value
#     db["elev"] = float("nan")
    db.loc[db["missing"]==True,"elev"] = None
    return db



def add_buffer(db,stats=None):
    #maybe rewrite the condition better
    if not(stats is None):
        if not np.isin(["buffer"],list(db)).any():
            raise Exception("Asked to add buffer to 'stats', 'db' must already have buffer")

        index = dplyr.filter(db, ~db["buffer"])
        index = index["seqno"]
        index["seqno_orig"] = list(range(len(index)))

        #TODO funny stuff at LITAP_load.R, add_buffer
        raise Exception('TODO funny stuff at LITAP_load.R, add_buffer')
        pass
        
    else:
        ncols= max(db["col"])
        nrows= max(db["row"])
        
        db.sort_values(by=["seqno"])
        
        if np.isin(["drec"],list(db)).any():
            drec = dplyr.distinct(db["drec"].to_frame()).reset_index(drop=True)
            
            
        names = ["row","col"]
        rows = [1] * (ncols+2) + list(range(1,nrows+3)) + list(range(1,nrows+3)) + [nrows+2] * (ncols+2)
        cols = list(range(1,ncols+3)) + [1] * (nrows+2) + [ncols+2] * (nrows+2) + list(range(1,ncols+3)) 
        buffers = pd.DataFrame(list(zip(rows,cols)),columns=names)
        buffers["buffer"] = True
        #maybe reset index after dropping duplicates? doing it further down and seems to be correct
        buffers = buffers.drop_duplicates()
        
        db["row"] = db["row"]+1
        db["col"] = db["col"]+1
        db["buffer"] = False
        
        
        db = pd.concat([db,buffers])
        db = db.sort_values(by=["row","col"])
        db = db.reset_index(drop=True)
        db["seqno_buffer"] = list(range(1,len(db["row"])+1))
        

        if np.isin(["drec"],list(db)).any():
            #TODO funny stuff at LITAP_load.R, add_buffer
            drec = dplyr.left_join(drec, dplyr.select(db,f.seqno,f.seqno_buffer),by = {'drec': 'seqno'})
            drec = drec.rename(columns={"seqno_buffer": "drec_buffer"})
            db = dplyr.left_join(db, drec,by = ['drec'])
            db = db.drop(columns=["drec"])
            db = db.rename(columns={"drec": "drec_buffer"})


            
            
        r = db.drop(columns=["seqno"])
        r = r.rename(columns={"seqno_buffer": "seqno"})
        names = r.columns.tolist()
        names = [names[-1]] + names[:-1]
        r = r[names]
        r = r.sort_values(by=["seqno"])
    return r


def db_prep(db, clim, rlim, edge, verbose):
    
    #if limits are provided
    if(rlim!=None and clim!=None):
        if(len(rlim)!=2 or len(clim)!=2):
            raise("clim and rlim must be each be a vector of two numbers (start/end row/col) or NULL")
        if((np.asarray(rlim)>max(db["row"])).any() or (np.asarray(clim)>max(db["col"])).any()):
            raise("Subset cannot be bigger than data")
        if((rlim[1]-rlim[0])<2 or (clim[1]-clim[0])<2):
            raise("Subset is too small (less than 2x2)")
        
        db = db[db["row"]>=rlim[0]] 
        db = db[db["row"]<=rlim[1]]
        db = db[db["col"]>=clim[0]] 
        db = db[db["col"]<=clim[1]]
        
        db["seqno"] = range(1, len(db) + 1)
        db["row"] = db["row"] - min(db["row"]) + 1
        db["col"] = db["col"] - min(db["col"]) + 1
       
    if(edge):
        
        db = add_buffer(db)
        
        #R code has once again sorting with seqno here, but it is alreay sorted from add_buffer function
        #backup of elevetions
        db["elev_orig"] =db["elev"]
        
        
        db1 = LITAP_functions.nb_values(copy.copy(db), max_cols = max(db["col"]), col = ["buffer", "missing"], format_ = "wide")
        
        
        buffer_cols = [x for x in db1.columns if x.startswith("buffer_")]
        missing_cols = [x for x in db1.columns if x.startswith("missing_")]
        
        
        db1["buffer_edge"] = db1.loc[:, db1.columns.isin(buffer_cols)].any(axis='columns')
        db1["missing_edge"] = db1.loc[:, db1.columns.isin(missing_cols)].any(axis='columns')
        
        db1["edge_map"] = db1.loc[:, db1.columns.isin(["buffer_edge","missing_edge"])].any(axis='columns')
        
        db = pd.merge(db,db1.loc[:, db1.columns.isin(["seqno","edge_map"])],how="left")
    return db


def load_file(file, nrow=None, ncol=None, missing_value=-9999,rlim=None,clim=None,edge=True,verbose=True):
    #TODO check if file exists
    
    #TODO get extension of file

    #TODO different loading for different files format
    
    db = load_dem(file,type_="elev")
    print()
    if (len(db.columns)>1 and np.isin(["x","y","elev"],list(db)).all()):
        #TODO preprocess if x,y is present as in R code (LITAP_load.R, laod_file)
        raise Exception('TODO in load_file')
    elif(nrow!=None and ncol!=None):
        print("using supplied nrow and ncol")
    else:
        raise Exception("dbf files with only one column require nrow and ncol arguments.")
        
    db = db_format(db, nrow, ncol, missing_value, verbose)
    
    db = db_prep(db, clim, rlim, edge, verbose)
    
    return db

def format_rule(rule,type_):
    for col in rule.columns:
        if rule[col].dtype == object:
            rule[col] = rule[col].str.lower()

    # !!!MAYBE THERE ARE MORE COLUMNS THAT HAVE slope IN THEIR NAME!!!
    if type_=="arule":
        rule["attr_in"] = rule["attr_in"].replace("slope", "slope_pct")
    
    if 'zone' not in rule.columns:
        rule["zone"] = 0

    return rule