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
from LITAP_functions import *
from LITAP_utils import *
from LITAP_load import *



#flow_calc_shed.R functions
def get_all_flow(db,w):
    
    db.loc[db["drec_shed"] == db["seqno_shed"],"drec_shed"] =  pd.NA
    m = db["seqno_shed"].to_frame()

    end = False
    t2 = 0
    
    while(not end):
        
        m1 = pd.Series(data=pd.NA, index=list(range(len(m))))
        
        last_column = m.iloc[:,-1]
        
        #find the positions of the fixed vector size in which to store the value
        valid_values_index = last_column[~pd.isna(last_column)].index
        
        #get the new indices who are stored as values
        m_last_col_indices = pd.Index(last_column[~pd.isna(last_column)].values-1)
        
        # update m1 vector with respect to the two stage indices mentioned before
        m1.iloc[valid_values_index] = db["drec_shed"].iloc[m_last_col_indices]
        
        indices_to_apply_NA =  m1[m1==last_column].index
        
        values_to_keep = pd.Series(data=pd.NA, index=list(range(len(m))))
        
        try:
            m1.iloc[indices_to_apply_NA] = pd.NA
        except Exception as e:
            print(t2)
            print(e)
            return values_to_keep, m1
        
        t2+=1
        m[str(t2)] = m1

        if m.iloc[:,-1].isnull().all():
            end = True

    return m

def get_upslope3(db, w, type_=["upslope","elev_diff"]):

    o = db["seqno_shed"]
#     db = db.sort_values(by=["seqno_shed"])#.reset_index(drop=True)  
    db = db >> dplyr.arrange(f.seqno_shed)
    
    
    if "elev_diff" in type_:
        db["elev_diff"] = 0
    if "upslope" in type_:
        db["upslope"] = 0

    
    if (w>0):

        flow = get_all_flow(copy.copy(db),w)

        t = [x for x in type_ if x in db.columns.tolist()]
        for cell in o:
            if ((db[t].iloc[cell-1]==0).any()):
                track = na_omit(flow.iloc[cell-1])

                if "upslope" in type_ and (db["upslope"].iloc[cell-1]==0):
                    db.loc[track-1,"upslope"] = upslope_values(track,db)
                if "elev_diff" in type_ and (db["elev_diff"].iloc[cell-1]==0):
#                     raise("TODO implement elev_diff_values function from flow_calc_shed.R")
                    db.loc[track-1,"elev_diff"] = elev_diff_values(track,db,w)
                    
    data_dict = {
        "shedno" : [w],
        "data" : [db]
    }

    return pd.DataFrame(data=data_dict)




def upslope_values(track,db):
    new = pd.Series(list(range(1,len(track)+1)))
    current = db.loc[track-1,"upslope"].reset_index(drop=True)

    try:
        new[current!=0] = new[current==0].iloc[-1]
    except:
        if new[current==0].empty:
            pass
        else:
            print("BAD")

    new = pd.Series(new.values,db.loc[track-1,"upslope"].index)
    
    return db.loc[track-1,"upslope"] + new


def elev_diff_values(track,db,w):
    new = dplyr.lag(db.loc[track-1, "elev"]) - db.loc[track-1, "elev"]
    new = new.drop(new.index[0])
    
    new = new.cumsum().reset_index(drop=True)
    
    current = db.loc[track-1, "elev_diff"]
    current = current.drop(current.index[-1]).reset_index(drop=True)
    
    try:
        new[current!=0] = new[current==0].iloc[-1]
    except:
        if new[current==0].empty:
            pass
        else:
            print("BAD")

    new = pd.Series(np.insert(new.values,0,0,axis=0),db.loc[track-1,"elev_diff"].index)            
    
    return db.loc[track-1,"elev_diff"] + new

def calc_upslopes(db, type_ = ["upslope","elev_diff"]):
    # Create sub-watershed seqno's for quicker referencing
    db = db.sort_values(by=["seqno"])#.reset_index(drop=True)   
    db["seqno_shed"] = db['shedno'].groupby(db['shedno'], group_keys=False).cumcount() + 1
    
    
    drec_index = db[~pd.isna(db["drec"])]["drec"].index
    drec_values_as_indices = pd.Index(db[~pd.isna(db["drec"])]["drec"].values-1)
    db["drec_shed"] = None
    db["drec_shed"] = db["drec_shed"].astype(pd.Int64Dtype())
    db["drec_shed"].loc[drec_index] = db["seqno_shed"][drec_values_as_indices].values
    
    # Calculate upslope flow for each cell by watershed
    db = db.sort_values(by=["elev"],ascending=False).reset_index(drop=True)  
    
    
    
    #calculate for each shedno in a parallel way the upslope
    db_cleared = db.replace(pd.NA,-9999)
    
    db_nested = db_cleared >> tidyr.nest(data=~f.shedno)
    db_nested = pd.concat(db_nested.apply(lambda x: get_upslope3(x["data"],x["shedno"],type_), axis=1).tolist())
    db_exploded = db_nested >> tidyr.unnest(f.data)
    db_exploded["drec_shed"] = db_exploded["drec_shed"].astype(pd.Int64Dtype())
    db_exploded = db_exploded.replace(-9999,pd.NA)
    
    db = db_exploded.sort_values(by=["seqno"]).reset_index(drop=True)
#     print(db.columns)
#     db = db.reindex(columns=["seqno","elev","drec","shedno"]+db.columns.tolist()[4:])
#     print(db.columns)
    
    
    return db


def calc_shed4(db, verbose = False):
    
    db_origin = copy.copy(db)
    
    npits = len(db["ddir"][db["ddir"]==5])
    
    db["shedno"] = None
    db["shedno"] = db["shedno"].astype(pd.Int64Dtype())
    
    indices = db[db["ddir"]==5].index
    numbers = list(range(1,npits+1))
    
    db["shedno"].loc[indices] = numbers
    
    db = db[["seqno","elev","drec","shedno"]]
    
    # Assign each cell to a watershed by climbing UP
    n1 = sum(~pd.isna(db["shedno"]))
    n2 = 0

    while(n1!=n2):
        drec_index = db[~pd.isna(db["drec"])]["drec"].index
        drec_values_as_indices = pd.Index(db[~pd.isna(db["drec"])]["drec"].values-1)
        db["shedno"].loc[drec_index] = db["shedno"][drec_values_as_indices].values
        
        n2 = n1
        
        n1 = sum(~pd.isna(db["shedno"]))
    
    # Relabel to match top down process
    
    shed_order = db.sort_values(by=["elev"],ascending=False).reset_index(drop=True)   
    shed_order = shed_order[~pd.isna(shed_order["shedno"])]
    shed_order = shed_order["shedno"]
    
#     db = 
    unique_pit_labels = shed_order.drop_duplicates().values
    initial_labels = list(range(1,npits+1))
    map_dict = dict(zip(unique_pit_labels,initial_labels))
    
    db = db.replace({"shedno": map_dict})

    # Calculate upslope and elev diff
    
    db = calc_upslopes(db,["upslope"])
#     print(db.iloc[4400:4403])
    # Save initial values
    db["initial_shed"] = db["shedno"]
    
    # Merge with original data
    db = db >> dplyr.left_join(db_origin,by=[f.seqno, f.elev, f.drec])
    db["drec"] = db["drec"].astype(pd.Int64Dtype())
    db["shedno"] = db["shedno"].astype(pd.Int64Dtype())
    db["initial_shed"] = db["initial_shed"].astype(pd.Int64Dtype())
    
    
    # Calculate ridge lines
    db1 = db >> dplyr.select(f.seqno, f.initial_shed, f.buffer)
    db1 = nb_values(db1, max_cols = max(db["col"]), col = ["initial_shed"])
    db1 = db1 >> dplyr.filter(f.buffer==False)

    db1 = db1.groupby(db1["seqno"]).apply(lambda x: len(np.unique(x.initial_shed_n[~pd.isna(x.initial_shed_n)])) > 1)
    db1_dict = {
        "seqno" : db1.index.tolist(),
        "ridge" : db1.values.tolist()
    }
    db1 = pd.DataFrame(db1_dict)
    db = db1 >> dplyr.right_join(db,by=[f.seqno])
    db = db >> dplyr.arrange(f.seqno)
    
    return db

