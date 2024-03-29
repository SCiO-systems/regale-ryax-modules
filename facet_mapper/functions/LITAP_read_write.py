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
import os

#import from other files
from LITAP_utils import *
from LITAP_load import *
import LITAP_functions   


def save_output(locs, out_format="csv", which=["fill","local", "pond", "pit", "ilocal"], where="flow",add_db=None):
    files = os.listdir(locs)
    
    for name in which:
        file_name = name + ".csv"
        if file_name in files:
            
            data = pd.read_csv(locs + file_name)

            if name in ["fill", "ilocal", "form", "weti", "relief", "length","fuzc", "fuza"] and "stats" not in name:
                
                if "db" in data.columns:
                    db = data["db"]
                else:
                    db = data
                
                if "data" in data.columns:
                    db = dplyr.select(db,~f.data)
                    
                if add_db!=None:
                    #TODO in save_output
                    raise("TODO in save_output")
                    
                
                save_shed(locs.replace("backup",where), db, "dem_" + file_name, clean=True)
            
            # I save stat files during the run of the script and not at the end from what I have saved in backup folder
            
            stat_names = ["fill","ilocal","initial","local","pit","pond"]
            if name in stat_names:
                #check if dataframe is empty
                stats_data = pd.read_csv(locs + "/stats_" + file_name)
                if len(stats_data>0):
                    s = remove_buffer(data,stats_data)
                    save_shed(locs.replace("backup",where),s,"stats_" + file_name)
            
            
    return 

def save_output2(data, locs, name, out_format, where, stats = None, add_db=None):
    
    if add_db is not None:
        data = dplyr.left_join(data, add_db, by = ["seqno"])
        
    if stats is not None:
        stats = remove_buffer(data,stats)
        save_shed(locs + where + "/", stats, "stats_" + name + ".csv")
    else:
        data = remove_buffer(data)
        save_shed(locs + where + "/", data, "dem_" + name + ".csv")
    
    return
    
def save_shed(file_out,obj,name,clean=False):
    
    if clean==True:
        # return obj
        obj = remove_buffer(obj)
        
    obj.to_csv(file_out + name,index=False)
    
    return

def remove_buffer(db, stats=None):
    index = dplyr.select(db,seqno_buffer= f.seqno)
    
    db = datar.all.filter(db,~f.buffer)
    db = datar.all.arrange(db,f.row,f.col)
    db = datar.all.rename(db, seqno_buffer = "seqno")
    
    if "drec" in db.columns:
        db = datar.all.rename(db, drec_buffer = "drec")
    
    row_len = len(db.row)
    db = dplyr.mutate(db,seqno=list(range(1,row_len+1)))
    
    for column_name in db.columns:
        if any(ext in column_name for ext in ["col","row"]):
            db[column_name] -=1
    
    index = dplyr.left_join(index,dplyr.select(db,f.seqno,f.seqno_buffer),by = ["seqno_buffer"])
    
    if stats is not None:
        #in R the command is: dplyr::mutate(dplyr::across(dplyr::matches("(^row)|(^col)|(_row)|(_col)"), ~ . - 1))
        for column in stats.columns:
            if "row" in column or "col" in column:
                stats[column] -= 1
                
        #in R the command is: dplyr::mutate(dplyr::across(dplyr::contains("seqno"), ~replace(., . == 0, as.numeric(NA)))) %>%
        for column in stats.columns:
            if "seqno" in column:
                stats[column] = stats[column].replace(0,pd.NA)

        #in R the command is: dplyr::mutate(dplyr::across(dplyr::contains("seqno"), ~index$seqno[.]))
        for column in stats.columns:
            if "seqno" in column:
                stats[column] = index.loc[stats[column]-1,"seqno"].reset_index(drop=True)
                
        return stats
    
    else:
        # Replace drec and upslope with correct cell numbers
        if "drec_buffer" in db.columns:
            db = dplyr.mutate(db, drec=index.loc[db["drec_buffer"]-1,"seqno"])
    return db

def get_previous(folder,step,where="backup",_type="dem"):
    #get folder items and try to find the wanted file 
    files_list = os.listdir(folder+where)
    files_list = [x for x in files_list if step in x]
    files_list = [x for x in files_list if _type in x]
    
    #check if conflicts, absense or all ok
    if len(files_list)>1:
        raise("more files than needed")
    elif len(files_list)==0:
        raise("not found the file")
    else:
        file_csv = pd.read_csv(folder + where + "/" + files_list[0])
        #if columns with _buffer exist drop them
        found_buffer_column = [x for x in file_csv.columns if "_buffer" in x]
        if found_buffer_column:
            file_csv = file_csv.drop(found_buffer_column,axis=1)
    return file_csv
