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
from flow_pit_stat_fast import *
from flow_calc_shed import *


# flow_remove_pit_fast.R functions
def remove_pit1(w_rm, w_stats, db, update_elev = False, verbose = False):
    
    w_rm = w_rm >> dplyr.mutate(direction = dplyr.if_else(f.pit_elev>=f.pit_elev_out,"out","in"))
    
    if w_rm["direction"][0]=="out":
        w_rm = w_rm >> dplyr.select(f.pit_seqno, f.in_seqno, f.out_seqno, f.pour_elev, f.shedno, f.drains_to)
    else:
        w_rm = w_rm >> dplyr.select(pit_seqno = f.pit_seqno_out, in_seqno = f.out_seqno, out_seqno = f.in_seqno,
                                    pour_elev = f.pour_elev_out, shedno = f.drains_to, drains_to = f.shedno)
    
    ## Remove Pit

    # Reverse directions of these cells (both regular and pit cell)

    # - Second reflows the pit cell. NOTE, if pit cell is the pour point, then
    #   this may be NA (fixed in next step)
    # - Third reflow the last cell in the new flow path (in_seqno) to point
    #   towards the first cell in the next flow path (out_seqno). This fixes any
    #   NA's from the previous step
    
        
        
    # 1. Get original flow down to old pit
    new_flow = trace_flow2(cell=w_rm["in_seqno"],drec=db["drec"])
    # Reverse directions
    new_ddir = db.iloc[new_flow] >> dplyr.mutate(ddir = 10 - dplyr.lag(f.ddir))
    # Fix pour_point direction to pour to new shed
    db_sub = new_ddir[new_ddir["seqno"]==w_rm["in_seqno"][0]].reset_index(drop=True)
    
    try:
        d = nb_values(db.drop(columns=["data"]), max_cols = max(db["col"]), col = ["seqno"], db_sub=db_sub.drop(columns=["data"]))
    except:
        d = nb_values(db, max_cols = max(db["col"]), col = ["seqno"], db_sub=db_sub)
        
        
    d = d >> dplyr.filter(f.seqno_n == w_rm["out_seqno"][0])
    d = d["n"][0]


    new_ddir.loc[new_ddir["seqno"]==w_rm["in_seqno"][0],"ddir"] = d
    new_ddir = flow_values(db, max_cols = max(db["col"]),col=["seqno"],db_sub=new_ddir)

     # Apply new flow directions

    db.ddir[new_ddir["seqno"]-1] = new_ddir["ddir"]
    db.drec[new_ddir["seqno"]-1] = new_ddir["seqno_next"]
    
    # Get new flow path from old pit centre to new pit centre
    new_flow = trace_flow2(cell=w_rm["pit_seqno"],drec=db["drec"])
    
    # Update upslope starting with old pit and looping through to new,
    #    adding upslope from previous to new cell along the way
    

    old_upslope = (db["upslope"].iloc[new_flow[:new_flow.index((w_rm["in_seqno"]-1)[0])+1]]).iloc[::-1].tolist()
    old_upslope.insert(0,0)

    lagged_list = (old_upslope- dplyr.lag(old_upslope,1)).iloc[::-1].tolist()
    del lagged_list[-1]
    lagged_list = [int(x) for x in lagged_list]
    
    new_upslope = np.cumsum(lagged_list).tolist()
    
    # Correct old shed upslope
    db["upslope"].iloc[new_flow[:new_flow.index((w_rm["in_seqno"]-1)[0])+1]] = new_upslope
    
    # Add to new shed upslope
    db["upslope"].iloc[new_flow[new_flow.index((w_rm["in_seqno"]-1)[0])+1:]] = db["upslope"].iloc[new_flow[new_flow.index((w_rm["in_seqno"]-1)[0])+1:]] + new_upslope[-1]
    
    # Update elevation of db (FlowMapR_2009.txt line 1964)
    if update_elev:
        indices = np.logical_and(
                    np.logical_and(db["shedno"]==w_rm["shedno"][0], db["elev"]<w_rm["pour_elev"][0]),
                    np.logical_and(~pd.isna(db["shedno"]), ~pd.isna(db["elev"])))
        db["elev"].loc[indices] = w_rm["pour_elev"][0]
    

    # Update watershed number
    db = db >> dplyr.mutate(shedno = dplyr.if_else(f.shedno==w_rm["shedno"][0],w_rm["drains_to"],f.shedno))
    
    return db

def first_pitr1(db, max_area=10,max_depth=0.5,verbose=False):
    
    # Working with initial_shed
    db["shedno"] = db["initial_shed"]
    
    # No point if nothing to remove
    if max_area>0 or max_depth>0:
        
        # Initial Shed Statistics
        w_stats = pit_stat1(db,verbose=verbose)
        w_stats = out_stat(w_stats)
        
        w_stats = w_stats >> dplyr.arrange(dplyr.desc(f.pit_elev),f.varatio)
        
        w_stats["remove"]= np.logical_and(np.logical_or(w_stats["pit_area"]<=max_area, (w_stats["pour_elev"]-w_stats["pit_elev"]<=max_depth)),~w_stats["edge_pit"])
        
        
    # In sequence, for each watershed that should be removed
    # Combine if too small (max_area) or too shallow (max_depth)
    done = False        
    while(not done):
        w_rm = w_stats >> dplyr.filter(f.remove==True)
        
        w_rm = w_rm >> dplyr.slice(1)
        if len(w_rm)>0:
            sheds = [w_rm["shedno"][0],w_rm["drains_to"][0]]
        
            db = remove_pit1(w_rm, sheds, db, update_elev=True)

            # Update shed statistics but only for the two sheds involved
            w_stats = w_stats >> dplyr.filter(~datar.base.testing.is_in(f.shedno,sheds))
            w_stats = w_stats >> dplyr.bind_rows(pit_stat1(db,w=sheds,verbose=verbose))
            w_stats = w_stats.dropna(subset=["shed_area"])
            # Which to remove
            w_stats["remove"]= np.logical_and(np.logical_or(w_stats["pit_area"]<=max_area, (w_stats["pour_elev"]-w_stats["pit_elev"]<=max_depth)),~w_stats["edge_pit"].astype(bool))
            # Sort order
            w_stats = w_stats >> dplyr.arrange(dplyr.desc(f.pit_elev),f.varatio)


            shed_rm = list(compress(sheds,~datar.base.testing.is_in(sheds,w_stats["shedno"])))[0] # Shed removed
            shed_kept = list(compress(sheds,datar.base.testing.is_in(sheds,w_stats["shedno"])))[0] # Shed removed            

            if any(~pd.isna(w_stats["drains_to"])):

                w_stats.loc[w_stats["drains_to"]==shed_rm,"drains_to"] = shed_kept
            
            w_stats = out_stat(w_stats)
                
                
        else:
            done = True
    # Update elev_diff
    db = calc_upslopes(db,type_=["elev_diff"])
    
    db["ridge"] = db["ridge"].replace(-9999,pd.NA)
    db["drec"] = db["drec"].replace(-9999,pd.NA)
    db["initial_shed"] = db["initial_shed"].replace(-9999,pd.NA)
    db["missing"] = db["missing"].replace(-9999,pd.NA)
    db["ddir"] = db["ddir"].replace(-9999,pd.NA)
    
    # Save as local_shed numbers
    db["local_shed"] = db["shedno"]
    db["local_ddir"] = db["ddir"]
    db["local_elev_diff"] = db["elev_diff"]
    

    return db


def second_pitr1(db, verbose = False):
    # Working with local_shed
    db["shedno"] = db["local_shed"]
    
    w_stats = pit_stat1(db,verbose=verbose)
    w_stats = out_stat(w_stats)
    w_stats = w_stats >> dplyr.arrange(f.pit_elev,f.varatio)
    
    w_stats = w_stats >> dplyr.mutate(removed=False, final = False, next_pit=f.drains_to, becomes=0)
    
    
    pond = w_stats
    removed = []
    final_pits = []
    done = False
    w = w_stats["shedno"].iloc[0]
    
    
    while(not done):
        
        w_rm = find_lowest(w, w_stats, final_pits=final_pits, removed=removed, verbose=verbose)
        
        w_focal = w_stats >> dplyr.filter(f.shedno == w_rm["shedno"].iloc[0])
        w_drain = w_stats >> dplyr.filter(f.shedno == w_rm["drains_to"].iloc[0])
        
        
        if not (w_focal["edge_pit"].iloc[0] and w_drain["edge_pit"].iloc[0]):
            
#             print(w_focal)
            new_shed = np.amax(db["shedno"]) + 1
#             print(new_shed)
#             print(db["shedno"])
#             new_shed = datar.base.arithmetic.max(db["shedno"],na_rm = True) + 1
            
            if verbose:
                print("  Combining sheds ", w_focal["shedno"].iloc[0], " and ", w_drain["shedno"].iloc[0], " into new shed ", new_shed)
        
            removed.append(w_rm["shedno"].iloc[0])
            removed.append(w_rm["drains_to"].iloc[0])
        
            
            # Remove shed
            db = remove_pit1(w_rm,w_stats,db)
        
        
            db.loc[datar.base.testing.is_in(db["shedno"].replace(pd.NA,-9999),[w_rm["shedno"].iloc[0],w_rm["drains_to"].iloc[0]]),"shedno"] = new_shed
            
                
            # Update shed statistics but only for the two sheds involved
            shed_update = datar.base.testing.is_in(w_stats["shedno"],[w_focal["shedno"].iloc[0],w_drain["shedno"].iloc[0]])

            w_stats.loc[shed_update,"removed"] = True
            w_stats.loc[shed_update,"final"] = w_rm["at_final"].iloc[0]
            w_stats.loc[shed_update,"next_pit"] = w_stats.loc[shed_update,"drains_to"]
            w_stats.loc[shed_update,"becomes"] = new_shed
            
            
            w_stats = w_stats >> dplyr.bind_rows(dplyr.mutate(
                pit_stat1(db,w=[new_shed],verbose=verbose),
                removed=False,
                final=False,
                next_pit=f.drains_to,
                becomes=0
            ))
            w_stats = w_stats.dropna(subset=['shed_area']).reset_index(drop=True)
            
            
            w_stats = w_stats >> dplyr.arrange(f.pit_elev,f.varatio)
            
            # Update references to old shed
            if(any( ~pd.isna(w_stats["drains_to"]))):
                w_stats.loc[datar.base.testing.is_in(w_stats["drains_to"],[w_focal["shedno"].iloc[0],w_drain["shedno"].iloc[0]]),"drains_to"] = new_shed
                
                w_stats.loc[np.logical_and(
                    np.logical_and(w_stats["removed"]==False,w_stats["final"]==False),
                    datar.base.testing.is_in(w_stats["next_pit"],[w_focal["shedno"].iloc[0],w_drain["shedno"].iloc[0]])),"next_pit"] = new_shed
            
            w_stats = out_stat(w_stats)
            
            # Keep track of pond statistics (ie. keep track of all watersheds)
            pond = dplyr.bind_rows(pond,w_stats.loc[w_stats["shedno"]==new_shed])
            
            sheds_update =  datar.base.testing.is_in(pond["shedno"],[w_focal["shedno"].iloc[0],w_drain["shedno"].iloc[0]])
            
            pond.loc[sheds_update,"removed"] = True
            pond.loc[sheds_update,"final"] = w_rm["at_final"].iloc[0]
            pond.loc[sheds_update,"next_pit"] = pond.loc[sheds_update,"drains_to"]
            pond.loc[sheds_update,"becomes"] = new_shed
            
            pond.loc[np.logical_and(
                np.logical_and(pond["removed"]==False,pond["final"]==False),
                np.logical_or(pond["next_pit"]==w_focal["shedno"].iloc[0],pond["next_pit"]==w_drain["shedno"].iloc[0])) , "next_pit"] = new_shed
            
            
            # Calc second vol2mm etc.
            
            
            
            vol = db >> dplyr.filter(f.shedno==new_shed)
            vol = vol >> dplyr.left_join(dplyr.select(w_stats,f.shedno,f.pour_elev,f.shed_area),by="shedno")
#             vol = vol >> dplyr.left_join(dplyr.select(w_stats,w_stats["shedno"],w_stats["pour_elev"],w_stats["shed_area"]),by="shedno")
            vol = vol2fl(vol,0, verbose=verbose)["data"].iloc[0]
            vol["shedno"] = new_shed
            
            # Only replace cells with new overflows (i.e. elev must be in vol)
            db_new = db >> dplyr.filter(f.shedno==new_shed, datar.base.testing.is_in(f.elev,vol["elev"]),f.parea==0)
            db_new = db_new >> dplyr.select(~f[f.vol2fl,f.mm2fl,f.parea])
            db_new = db_new >> dplyr.left_join(vol,by = ["shedno","elev"])
            db_new = db_new >> dplyr.arrange(f.seqno)
            
            
            db_new.loc[ ~pd.isna(db_new["parea"]), ["mm2fl","vol2fl","parea"]] = 0
            
#             print(db_new["vol2fl"])
            
#             print(db.loc[db_new["seqno"]-1, ["vol2fl","parea"]])
#             print( db_new.loc["vol2fl","parea"])
            
            db.loc[db_new["seqno"]-1, ["vol2fl","parea"]] = db_new[["vol2fl","parea"]].values
            
        
        else:
            
            if(verbose):
                print("  Watersheds ", w_focal["shedno"].iloc[0], " and ", w_drain["shedno"].iloc[0], " are FINAL sheds")

            final_pits =  list(np.unique(final_pits+w_rm["shedno"].tolist()+w_rm["drains_to"].tolist()))
            
            w_stats.loc[datar.base.testing.is_in(w_stats["shedno"],[w_focal["shedno"],w_drain["shedno"]]) , "final"] = True
            
            
            
        # Get next shed
        finished = list(np.unique(final_pits+removed))

        w = w_stats >> dplyr.filter( ~datar.base.testing.is_in(f.shedno,finished))

        if len(w)>0:
            w = w["shedno"].iloc[0]
        else:
            done = True
    
    # Update final sheds
    
    pond.loc[datar.base.testing.is_in(pond["shedno"],final_pits),"final"] = True
    
    # Save as pond_shed numbers
    db["pond_shed"] = db["shedno"]
    
    db = {"db" : db,
          "stats" : pond}
    
    return db

def third_pitr1(db, verbose = False):
    
    # Working with local_shed
    db["shedno"] = db["local_shed"]

    w_stats = pit_stat1(db,verbose=verbose)
    w_stats = out_stat(w_stats)
    w_stats = w_stats >> dplyr.arrange(f.varatio)
    
    w_stats = w_stats >> dplyr.mutate(drains_to_orig = f.drains_to, next_pit = f.drains_to, end_pit=f.pond_shed, removed=False)
    
    fill = pd.DataFrame()
    done = False
    finished = []
    w = w_stats["shedno"].iloc[0]
    
    
    while(not done):
        w_focal = w_stats >> dplyr.filter(f.shedno==w)
        
        w_drain = w_stats >> dplyr.filter(f.shedno==w_focal["drains_to"].iloc[0])
        
        if w_focal["end_pit"].iloc[0]==w_drain["end_pit"].iloc[0]:
            new_shed = np.amax(db["shedno"]) + 1
            
            if verbose:
                print("  Combining sheds ", w_focal["shedno"].iloc[0], " and ", w_drain["shedno"].iloc[0], " to new shed ", new_shed)
                

            # Remove shed
            db = remove_pit1(w_focal,w_stats,db)
        
        
            db.loc[datar.base.testing.is_in(db["shedno"].replace(pd.NA,-9999),[w_focal["shedno"].iloc[0],w_drain["shedno"].iloc[0]]),"shedno"] = new_shed
            
            # Update shed statistics but only for the two sheds involved
            
            w_new = pit_stat1(db, w = [new_shed], verbose=verbose)
            w_new = w_new >> dplyr.mutate(drains_to_orig = f.drains_to, next_pit = f.drains_to, end_pit = f.pond_shed, removed = False)

            #check if ~ needs not instead
            w_stats = w_stats >> dplyr.filter(~datar.base.testing.is_in(w_stats["shedno"].replace(pd.NA,-9999),[w_focal["shedno"].iloc[0],w_drain["shedno"].iloc[0]]))
            
            w_stats = w_stats >> dplyr.bind_rows(w_new)
            w_stats = out_stat(w_stats)
            w_stats = w_stats.dropna(subset=['shed_area']).reset_index(drop=True)
            # Sort order
            w_stats = w_stats >> dplyr.arrange(f.varatio)
            
            # Update references to old shed
            if(any( ~pd.isna(w_stats["drains_to"]))):
                w_stats = w_stats >> dplyr.mutate(drains_to_orig = f.drains_to)
                
                w_stats.loc[datar.base.testing.is_in(w_stats["drains_to"],[w_focal["shedno"].iloc[0],w_drain["shedno"].iloc[0]]),"drains_to"] = new_shed
            
            w_stats = out_stat(w_stats)
            
            
            # Update fill data frame
            # check if this is correct
            w_focal["becomes"] = new_shed
            w_drain["becomes"] = new_shed
            temp = dplyr.bind_rows(w_focal,w_drain)
            temp["final"] = False
            temp = temp >> dplyr.arrange(f.shedno)
            fill = fill >> dplyr.bind_rows(temp)
            
            
            finished.append(w_focal["shedno"].iloc[0])
            finished.append(w_drain["shedno"].iloc[0])
            
            
            # Update mm2fl
            db_new = db >> dplyr.filter(f.shedno==new_shed, f.vol2fl!=0)
            db_new = db_new >> dplyr.left_join(dplyr.select(w_stats,f.shedno,f.shed_area),by=f.shedno)
            db_new = db_new >> dplyr.filter(f.mm2fl < w_focal["varatio"].iloc[0])
            db_new = db_new >> dplyr.mutate(mm2fl = f.vol2fl/f.shed_area)
            
            
            db.loc[db_new.loc[~pd.isna(db_new["seqno"]),"seqno"] - 1,"mm2fl"] = db_new.loc[~pd.isna(db_new["seqno"]),"mm2fl"]
        
        else:
            finished.append(w_focal["shedno"].iloc[0])
            
        
        # Get next shed
        w = w_stats >> dplyr.filter( ~datar.base.testing.is_in(f.shedno,finished))

        if len(w)>0:
            w = w["shedno"].iloc[0]
        else:
            done = True
        
        
    # When all is done, add sheds that weren't removed
    add = w_stats

    if len(fill)>0:

        add = add >> dplyr.filter(~datar.base.testing.is_in(f.shedno,fill["shedno"]) and f.removed==False)

    add["becomes"] = add["shedno"]
    add["final"] = True
    add = add >> dplyr.arrange(f.shedno)


    fill = fill >> dplyr.bind_rows(add)
    fill["next_pit"] = fill["drains_to"]
    fill["drains_to"] = fill["drains_to_orig"]
    fill = fill >> dplyr.select(~f.drains_to_orig)
    
    db["fill_shed"] = db["shedno"]
    
    # Update elev_diff
    if verbose:
        print("    Calculating new elevation differences")
        
    db = calc_upslopes(db,type_=["elev_diff"])
    
    
    db["ridge"] = db["ridge"].replace(-9999,pd.NA)
    db["drec"] = db["drec"].replace(-9999,pd.NA)
    db["initial_shed"] = db["initial_shed"].replace(-9999,pd.NA)
    db["missing"] = db["missing"].replace(-9999,pd.NA)
    db["ddir"] = db["ddir"].replace(-9999,pd.NA)
    
    db["elev"] = db["elev"].replace(-9999,pd.NA)
    db["elev_orig"] = db["elev_orig"].replace(-9999,pd.NA)
    db["local_shed"] = db["local_shed"].replace(-9999,pd.NA)
    db["local_ddir"] = db["local_ddir"].replace(-9999,pd.NA)
    db["local_elev_diff"] = db["local_elev_diff"].replace(-9999,pd.NA)
    db["pond_shed"] = db["pond_shed"].replace(-9999,pd.NA)
    db["fill_shed"] = db["fill_shed"].replace(-9999,pd.NA)
    
    db = {"db" : db,
      "stats" : fill}
    
    return db
    
