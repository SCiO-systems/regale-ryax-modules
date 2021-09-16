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


#flow_pit_stat_fast.R functions
def pit_stat1(db, w = None, verbose = False):
    
    # If more than one watershed
    if (len(np.unique(db.shedno[~pd.isna(db.shedno)])) > 1):
        
        # Subset watersheds
        if w is None:
#             print(db["shedno"])
            w = np.unique(db.shedno[~np.isnan(db["shedno"].astype(pd.Int64Dtype()))])
        
        
        # For each watershed calculate the pour_point (details of the point at which tip to another watershed)
        
        # Take only ridge cells
        db_pits = db[np.logical_and(~pd.isnull(db.shedno),db["ridge"]==True)]
        db_pits = db_pits[["seqno", "shedno", "elev", "upslope"]].reset_index(drop=True)
    
        db_pits = nb_values(db, max_cols = max(db["col"]), col =["elev", "shedno", "seqno"],db_sub=db_pits)
        
        
        pp = db_pits >> dplyr.filter(f.shedno!=f.shedno_n)
        pp = pp.dropna(subset=['shedno_n']).reset_index(drop=True)
        pp["seqno_n"] = pp["seqno_n"].astype(pd.Int64Dtype())
        pp["shedno_n"] = pp["shedno_n"].astype(pd.Int64Dtype())
        
        
        if (len(np.unique(pp["shedno"]))!=len(np.unique(db_pits.shedno[~pd.isnull(db_pits.shedno)]))):
            raise("Is implemented but need to be checked with data, flow_pit_stat_fast.R")
            pp_temp = db_pits[np.logical_and(~(db_pits["shedno"].isin(pp["shedno"])),db_pits["shedno"]!=db_pits["shedno_n"])]
            pp = pp_temp >> dplyr.bind_rows(pp)
        
        
        pp["pour_elev"] = pp[["elev","elev_n"]].apply(lambda x: max(x), axis=1)
        
        pp["min_elev_n"] = pp["elev_n"].groupby(pp["seqno"]).transform(lambda x: min(x))
        
        pp = pp[pp["elev_n"]==pp["min_elev_n"]]
        
        pp = pp >> dplyr.group_by(f.shedno)
        pp = pp >> dplyr.filter(f.pour_elev==datar.all.min(f.pour_elev))
#         pp = pp >> dplyr.ungroup()
        
        pp = pp >> dplyr.arrange(f.shedno,f.seqno)
        pp = pp >> dplyr.slice(1)
        
        
        
        pp = pp >> dplyr.left_join(dplyr.select(db,f.seqno,f.col,f.row),by = ["seqno"])
        
        pp["col_out"] = pp["seqno_n"].apply(lambda x: db["col"].iloc[x-1])
        pp["row_out"] = pp["seqno_n"].apply(lambda x: db["row"].iloc[x-1])
        
        pp = pp >> dplyr.select(f.shedno, f.pour_elev, in_seqno = f.seqno, in_row = f.row, in_col = f.col, in_elev = f.elev,
                    out_seqno = f.seqno_n, out_row = f.row_out, out_col = f.col_out, out_elev = f.elev_n, drains_to = f.shedno_n)
        
#         print(db)
        # Add other calculations
        stats = db[db["shedno"].isin(w)]
#         print(db["shedno"].isin(w))
#         print(pd.Series(datar.base.testing.is_in(db["shedno"],w)))
#         stats = db.loc[datar.base.testing.is_in(db["shedno"],w)]
#         print(stats)
        stats = stats >> dplyr.left_join(pp,by=f.shedno)
        stats = stats >> dplyr.group_by(f.shedno)
        stats = stats >> dplyr.mutate(shed_area = datar.base.length(f.shedno))
        stats = stats >> dplyr.filter(f.elev<=f.pour_elev)
        
        stats = stats >> dplyr.summarise(
                                shed_area = f.shed_area[0],
                                edge_pit = datar.base.testing.any(f.edge_map),
                                pit_area = datar.base.length(f.shedno),
                                pit_vol = datar.base.arithmetic.sum(f.pour_elev - f.elev),
                                pit_elev = f.elev[f.ddir == 5],
                                pit_seqno = f.seqno[f.ddir == 5],
                                pit_row = f.row[f.ddir == 5],
                                pit_col = f.col[f.ddir == 5],
                                pre_vol = 0,
                                varatio = dplyr.if_else(f.shed_area > 0,
                                f.pit_vol / f.shed_area * 1000, 0)) 
        
        stats = stats >> dplyr.right_join(pp,by=f.shedno)
    
    else:
        raise("TODO else bracket in flow_pit_stat_fast.R")
    

    if ((db.columns).isin(["pond_shed"]).any()):
        
        temp = db >> dplyr.select(f.shedno,f.pond_shed)
        temp = temp >> dplyr.distinct()
        stats = stats >> dplyr.left_join(temp, by=f.shedno)
        
#         raise("TODO small if bracket in flow_pit_stat_fast.R")
    
    return stats



def out_stat(pit_stats):
    
    
    
    pit_stats = pit_stats >> dplyr.select(~dplyr.ends_with("out"))
#     pit_stats["shedno_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats["shedno"].iloc[x-1])
#     pit_stats["edge_pit_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats["edge_pit"].iloc[x-1])
    
#     print(pit_stats)
    pit_stats["edge_pit_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats.loc[pit_stats["shedno"]==x,"edge_pit"]).bfill(axis=1).iloc[:,0]
#     print(pit_stats["drains_to"].apply(lambda x: pit_stats.loc[pit_stats["shedno"]==x,"edge_pit"]).bfill(axis=1)[[0]])

#     pit_stats["pit_elev_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats["pit_elev"].iloc[x-1])
    pit_stats["pit_elev_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats.loc[pit_stats["shedno"]==x,"pit_elev"]).bfill(axis=1).iloc[:,0]
    
#     pit_stats["pit_seqno_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats["pit_seqno"].iloc[x-1])
    pit_stats["pit_seqno_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats.loc[pit_stats["shedno"]==x,"pit_seqno"]).bfill(axis=1).iloc[:,0]
    
#     pit_stats["pour_elev_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats["pour_elev"].iloc[x-1])
    pit_stats["pour_elev_out"] = pit_stats["drains_to"].apply(lambda x: pit_stats.loc[pit_stats["shedno"]==x,"pour_elev"]).bfill(axis=1).iloc[:,0]
    
    
    return pit_stats



def vol2fl(db, w, verbose):
    if (db["shed_area"]<=0).any():
        raise("Shed area <= 0, is this reasonable?")
    
    vol_stats = db >> dplyr.arrange(f.elev,dplyr.desc(f.upslope))
    vol_stats = vol_stats >> dplyr.filter(f.elev<=f.pour_elev)
    vol_stats = vol_stats >> dplyr.group_by(f.elev,f.shed_area)
    vol_stats = vol_stats >> dplyr.summarise(total_cells = datar.all.length(f.elev))
    vol_stats = vol_stats >> dplyr.ungroup()
    vol_stats = vol_stats >> dplyr.mutate(parea = datar.all.cumsum(f.total_cells),
                                          last_elev = dplyr.lag(f.elev, default = f.elev[0]),
                                          elev_diff = (f.elev - f.last_elev) * 1000,
                                          vol2fl = pd.NA)
    
        
        
        
    if len(vol_stats)>0:
        prev = pd.NA
        for i in range(0,len(vol_stats)):
        
            if prev is pd.NA:
                prev = 0.1
                
            vol_stats["vol2fl"].iloc[i] = prev + vol_stats["elev_diff"].iloc[i] * (vol_stats["parea"].iloc[i]-1)
            
            prev = vol_stats["vol2fl"].iloc[i]
        
        
    
    vol_stats = vol_stats >> dplyr.mutate(mm2fl = dplyr.if_else(f.shed_area>0,f.vol2fl/f.shed_area,f.vol2fl/1))
    vol_stats = vol_stats >> dplyr.select(f.elev,f.vol2fl,f.mm2fl,f.parea)
    
    
    data_dict = {
        "shedno" : [w],
        "data" : [vol_stats]
    }

    return pd.DataFrame(data=data_dict)  


def calc_vol2fl(db, i_stats, verbose):
    # According to DEMProces.cpp (C++ flowmapr), Edge pits are skipped
    
    i_stats = i_stats >> dplyr.filter(f.edge_pit==False)
    
    db = db >> dplyr.mutate(shedno = f.initial_shed)
    
    # at R source code there is a command that gets rid of any existing R lists in the dataframe
    # dont know if there is a nessecity for this command in python, will skip it for the moment
    
    vol = db >> dplyr.filter(~datar.base.na.is_na(f.shedno))
    
    if (len(vol)>0 and len(i_stats)>0):
        
        vol = vol >> dplyr.right_join(dplyr.select(i_stats,f.shedno,f.pour_elev,f.shed_area),by=["shedno"])
        vol = vol >> tidyr.nest(data=~f.shedno)
        vol["vol"] = pd.concat(vol.apply(lambda x: vol2fl(x["data"],x["shedno"],verbose), axis=1).tolist()).reset_index(drop=True).data
#         return vol
        vol = vol >> tidyr.unnest(f.vol)
        
        
        db["shedno"] = db["shedno"].replace(pd.NA,-9999)
        
        db = dplyr.left_join(db,vol,by=[f.shedno,f.elev])
        
        db["shedno"] = db["shedno"].replace(-9999,pd.NA)
        
        db.loc[datar.base.na.is_na(db["parea"]),["mm2fl","vol2fl","parea"]] = 0
        
    else:
        db = db >> dplyr.mutate(vol2fl=0,mm2fl=0,parea=0)
        
    return db