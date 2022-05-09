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
from LITAP_utils import *
from LITAP_load import *


#LITAP_functions.R functions 
def nb_values(db, max_cols, col = ["elev"], db_sub = None, format_ = "long"):
    
    if db_sub is None:
        db_sub = db
        
    for i in range(1,10):
        if (i==1):
            seqno = db_sub["seqno"] + (max_cols - 1)
        elif (i==2):
            seqno = db_sub["seqno"] + max_cols
        elif (i==3):
            seqno = db_sub["seqno"] + (max_cols + 1)
        elif (i==4):
            seqno = db_sub["seqno"] - 1
        elif (i==5):
            seqno = db_sub["seqno"] 
        elif (i==6):
            seqno = db_sub["seqno"] + 1
        elif (i==7):
            seqno = db_sub["seqno"] - (max_cols + 1)
        elif (i==8):
            seqno = db_sub["seqno"] - max_cols
        elif (i==9):
            seqno = db_sub["seqno"] - (max_cols - 1)
    
        #copy.copy for the usual problem of pointers in python in multistage structures
        seqno =copy.copy(seqno) - 1
              
        max_seqno = max(seqno)
        min_seqno = min(seqno)
        
        max_db_seqno = max(db["seqno"])-1
        min_db_seqno = min(db["seqno"])-1
        
        keep_list = np.logical_and(seqno>=0,seqno<=max_db_seqno).to_list()
        seqno_list = seqno.iloc[keep_list].tolist()

  

        for a in col:
            if (min_db_seqno-min_seqno)>0:
                r = db[a].iloc[seqno_list].reindex(range(min_seqno-min_db_seqno,len(seqno_list))).reset_index(drop=True)

            elif(max_db_seqno-max_seqno<0):
                r = db[a].iloc[seqno_list].reset_index(drop=True).reindex(range(0,len(seqno_list)+abs(max_db_seqno-max_seqno)))

            else:
                #using t to keep the dimensions right. Problem showed up at flatcells because their total number was smaller than the length of the db
                t = np.zeros(len(db))
                t[seqno_list] = 1
                t = t> 0
                r =  db[a].iloc[t].reset_index(drop=True)
            
            
            db_sub[a+"_n"+str(i)] = r

    if format_=="long": 
        reg_expr = "("
        for i,column_name in enumerate(col):
            if i==0:
                reg_expr = reg_expr + column_name + "_n[0-9]{1})"
            else:
                reg_expr = reg_expr + "|(" + column_name + "_n[0-9]{1})"
        
        stubnames = []
        for matches in db_sub.columns.str.findall(reg_expr).values:
            if matches == []: 
                continue
            else:
                if (type(matches[0])!=str):
#                     print("diaforo listas")
                    matches = matches[0]

                for match in matches:
                    if match=="":
                        continue
                    else:
                        stubnames.append(match)
        
        ##tiring procedure for making wide to long
        # creating an id column with the index values
        db_sub["id"] = db_sub.index
        #changing from wide to long
        db_sub = pd.melt(db_sub,id_vars=[item for item in db_sub.columns if item not in stubnames],value_vars=stubnames,var_name="n")
        #splitting one column to two and doing some "make-up" to them
        db_sub[['type', 'n']] = db_sub['n'].str.split('_n', 1, expand=True)
        db_sub["type"] = db_sub["type"] + '_n'
        db_sub["n"] = (db_sub["n"]).astype(int)
        #substituting NaN values to -9999 in order for the pandas funtions to work
        db_sub = db_sub.fillna(-9999)
        #spreading the dataframe based on the values of one column
        db_sub = db_sub.pivot_table(index = [item for item in db_sub.columns if item not in ["value","type"]], values = 'value', columns = 'type').reset_index()
        
        # have in mind if the sorting needs changing, inlcuding the "col"
        #doing some "make-up" on the dataframe to bring it to the correct form
        #sort firstly with the neighbour n and then with the global id in order to get the right order in the final dataframe
        db_sub = db_sub.sort_values(by=["n","id"]).reset_index(drop=True)
        db_sub.columns.name = None
        db_sub = db_sub.replace(-9999.0,np.nan)
        db_sub = db_sub.drop(["id"], axis=1)

    
    return db_sub
    
    
def finddir2(db):
    #find difference in elevation between central pixel and all its neighbors
    db["elev_diff"] = db["elev"] - db["elev_n"]
    
    #adjust the distance to the edge pixels with sqrt(2)
    db.loc[db["n"]%2==1,'elev_diff'] = db.loc[db["n"]%2==1,'elev_diff']/np.sqrt(2)
    
    #replace nan values with numerical values of 0, 
    db["elev_diff"] = db["elev_diff"].groupby(db['seqno'], group_keys=False).apply(lambda x: x.replace(np.nan,0))
    #find the maximum elevation difference of the central pixels with its neighbors
    db["max_slope"] = db["elev_diff"].groupby(db['seqno'], group_keys=False).transform('max')
    
    
    
    #mark the direction of the flow towards the maximum difference neighbor
    #it is done in 2 steps because R-Python compatibility problems or lack of knowledge and motivation to find a better solution
    #in the first step the direction is calculated
    db["ddir"] = db.groupby(db['seqno'], group_keys=False).apply(lambda x: x.loc[np.logical_and(x["elev_diff"]==x["max_slope"], x["max_slope"]>0),"n"])
    
    #in the second step this value is written in the ddir column of all the neighbors of the central pixel
    def my_func(t):
        if t[~pd.isna(t)].empty:
            return pd.NA
        else:
            return t[~pd.isna(t)].iloc[0]
        
    db["ddir"] = db["ddir"].groupby(db['seqno'], group_keys=False).transform(lambda x: my_func(x))

    db = db.drop(["elev_n","max_slope","n","elev_diff"], axis=1)
    
    db = db.drop_duplicates()
    
    db.loc[np.logical_and(pd.isna(db["ddir"]), ~pd.isna(db["elev"])),"ddir"] = 5
    
    db["flatcell"] = db["ddir"]==5
    
    return db

def flow_values(db, max_cols, col = ["elev"], db_sub = None):

    if db_sub is None:
        db_sub = db

    for a in col:
        f = db_sub["seqno"]
        g = db_sub["ddir"].apply(lambda x: calc_seq(x,max_cols))
#         print(f.iloc[4398:4404])
#         print(g.iloc[4398:4404])
#         print(db_sub["ddir"].iloc[4398:4404])
        indexing = f + g
        db_sub[a + "_next"] = np.nan
        t = indexing[~pd.isna(indexing)].astype(int) - 1
        db_sub[a + "_next"].loc[t.index] = db[a][t.tolist()].values
        db_sub["seqno_next"] = db_sub["seqno_next"].astype(pd.Int64Dtype())
        db_sub["ddir"] = db_sub["ddir"].astype(pd.Int64Dtype())
    return db_sub

def calc_seq(d, max_cols):
    if d is pd.NA:
        return pd.NA
    if (d==1):
        return max_cols - 1
    elif (d==2):
        return max_cols
    elif (d==3):
        return max_cols + 1
    elif (d==4):
        return -1
    elif (d==5):
        return 0
    elif (d==6):
        return 1
    elif (d==7):
        return -(max_cols + 1)
    elif (d==8):
        return -max_cols
    elif (d==9):
        return -(max_cols - 1)

def trace_flow2(cell, drec):
#     print(drec.iloc[4400])
#     print("******************************")
    
    track = []
    if len(cell)>0:
        cell = int(cell[0]) - 1
        track = [cell]
        end = False
        while(not end):
#             print(cell)
#             print(drec.iloc[cell+1])
            cell = drec.iloc[cell] -1
#             print(cell)
#             print("----")
            if (not (cell is pd.NA)):
                if cell in track:
#                     print("cell:",cell," is in list")
                    end = True
                if cell!=track[-1]:
                    track.append(cell)
            else:
                end = True
#     print([x+1 for x in track])
    return track


def trace_pits(shedno, w_stats):
    
    track = [shedno]
    end =  False
    
    while(not end):
        shedno = w_stats.loc[w_stats["shedno"]==shedno,"drains_to"].iloc[0] # Get next cells
        
        
        if datar.base.testing.is_in(shedno, track): 
            end = True # In a circular track
        if shedno != track[-1]:
            track.append(shedno) # If not simply starting and ending with the same shed
    return track


def find_lowest(w, w_stats, final_pits, removed, verbose=False):
    # # Mark final pits
    # w_stats = dplyr::mutate(w_stats,
    #                          final = shedno %in% final_pits,
    #                          removed = shedno %in% removed)

    # Starting point
    
    w_pits = trace_pits(w, w_stats)
    
    wp = w_pits[0]
    lowest = w_stats >> dplyr.filter(f.shedno==wp)
    
    visited = []
    end = False
    
    while(not end):

        if(not datar.base.testing.is_in(wp,visited)):
            visited.append(wp)
            
            pit1 = w_stats >> dplyr.filter(f.shedno==wp)
            pit1["at_final"] = False
        

            pit2 = w_stats >> dplyr.filter(f.shedno==pit1["drains_to"].values[0])
            pit2["at_final"] = False
            
        
        
            
            if pit1["pit_elev"].iloc[0] < lowest["pit_elev"].iloc[0]:
                lowest = pit1
                
                
            if verbose:
                print("  Current pit: ", wp)
                print("    Pit 1: ", pit1["shedno"].iloc[0])
                print("    Pit 2: ", pit2["shedno"].iloc[0])
            
#             print(pit2)
            wp = pit2["shedno"]
            
            if pit2["final"].iloc[0]:
                lowest = pit1
                end = True
                if verbose:
                    print("    Pit 2 already FINAL pit")
            else:
                if datar.base.testing.is_in(pit2["shedno"].iloc[0],visited):
                    end = True
                else:
                # Technically this line here makes it impossible to go more than one pit down...
                    visited.append(pit2["shedno"])
                    
                    if pit2["pit_elev"].iloc[0] < lowest["pit_elev"].iloc[0]:
                        lowest = pit2
                        
                    if pit2["pour_elev"].iloc[0] < pit1["pour_elev"].iloc[0]:
                        lowest = pit2
                    else:
                        if lowest["final"].iloc[0]:
                            lowest = w_stats >> dplyr.filter(f.shedno == w_pits[0])
                            lowest["at_final"] = True
                        else:
                            lowest["at_final"] = False
            
                        end = True
                
        else:
            end = True
        
    
    if verbose:
        print("    Lowest pit: ", lowest["shedno"].iloc[0])
    
    return lowest



def get_dir(row, col, row_f, col_f, ddir_opts = list(range(1,10))):
#     print(row, col, row_f, col_f, ddir_opts)
    h = col_f-col
    v = row_f -row
    if v!=0:
        a = math.atan(abs(h)/abs(v)) * 180/math.pi
    else:
        a = math.pi/2 * 180/math.pi
        
    
    if h==0 and v==0:
        l = 5
    elif h>0 and v<=0: # Move to Upper right
        if a>=67.5:
            l = 6
        elif a>=22.5:
            l = 9
        else:
            l = 8
    elif h>0 and v >0: # Lower right
        if a>=67.5:
            l = 6
        elif a>=22.5:
            l = 3
        else:
            l = 2
    elif h<=0 and v <=0: # Upper left
        if a>=67.5:
            l = 4
        elif a>=22.5:
            l = 7
        else:
            l = 8
    elif h<=0 and v>=0: # Lower left
        if a>=67.5:
            l = 4
        elif a>=22.5:
            l = 1
        else:
            l = 2            
            
            
  # If best direction not an option, get next best
  # Prioritorize smaller seqno
    if not datar.base.testing.is_in(l,ddir_opts) and l!=5:
        raise("TODO: if case in get_dir function in LITAP_functions.R, called from calc_dir2 function in flow_calc_dir.R file. Relatively easy stuff.")
        closest = {
            "1" : [4, 2, 7, 3, 8, 6, 9],
            "2" : [1, 3, 4, 6, 7, 9, 8],
            "3" : [6, 2, 9, 1, 8, 4, 7],
            "4" : [7, 1, 8, 2, 9, 3, 6],
            "6" : [9, 3, 8, 2, 7, 1, 4],
            "7" : [8, 4, 9, 1, 6, 2, 3],
            "8" : [7, 9, 4, 6, 1, 3, 2],
            "9" : [8, 6, 7, 3, 4, 2, 1]
        }
    
        c = closest[str(l)][0]
    
    return l 

def trace_flow_all(cells, drec):
    m = np.expand_dims(cells.values,axis=0)
    i = 0
    while True:
        i += 1
        if all(m[i-1,:]==drec.iloc[m[i-1,:]-1].values.astype(int)):
            break
        m = np.concatenate((m,np.expand_dims(drec.iloc[m[i-1,:]-1].values,axis=0)),axis=0)
    
    return m
