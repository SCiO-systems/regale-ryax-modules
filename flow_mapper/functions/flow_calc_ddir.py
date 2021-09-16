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


#flow_calc_ddir.R functions
def calc_ddir2(db,verbose=False):
    # Calculate elevation of 'shifted' directions (neighbours)
    db = nb_values(db,max_cols = max(db["col"]),col = ["elev"])
    
    db1 = finddir2(db) # Calc Flow direction
    # Check for flat plateaus
    # - get neighbouring directions

    db_flats = db1 >> dplyr.filter(f.flatcell == True)
    
    db_flats = nb_values(db1, max_cols=max(db["col"]), col = ["ddir","elev"], db_sub = db_flats.reset_index(drop=True))
    db_flats = db_flats >> dplyr.group_by(f.seqno)
    
    end = False
    a = 0
    
    
    #Much things to be done here 
    while(not end):
        
        # Calculate change in ddir
        db_flats = db_flats >> dplyr.mutate(elev_diff = f.elev[0]-f.elev_n)
        
        db_flats = db_flats >> dplyr.summarise(ddir = f.n[(f.elev_diff>=0) & (f.elev_diff!=pd.NA) & (f.ddir_n!=5) & (f.ddir_n!=pd.NA)])
        db_flats = db_flats >> dplyr.ungroup()
        db_flats = db_flats.dropna().drop_duplicates(subset=["seqno"])

        

        if 'ddir' not in db_flats.columns:
            end = True
            
        if not end:
            db_flats = db_flats.loc[~pd.isna(db_flats['ddir'])]
            if len(db_flats)==0:
                end = True
            else:
                # Change db1
#                 print(db_flats)
                db1.loc[(db_flats["seqno"]-1).tolist(),"ddir"] = db_flats["ddir"].values
                
                # Recalculate neighbours
                db_flats = db1.loc[np.logical_and(~pd.isna(db1["ddir"]),db1["ddir"]==5)].reset_index(drop=True)
                db_flats = nb_values(db1, max_cols=max(db["col"]), col = ["ddir","elev"], db_sub = db_flats.reset_index(drop=True))
                db_flats = db_flats >> dplyr.group_by(f.seqno)
        

    
          
    
    # Deal with flow into depressions (flatin)
    # - Get the area of the depression, assign middle cell to a slightly lower elevation to break ties
    
    # Get patches of pit cells    
#     db_flats = db1.loc[np.logical_and(~pd.isna(db1["ddir"]),db1["ddir"]==5)].reset_index(drop=True)
    db_flats = db1 >> dplyr.filter(f.ddir==5,f.ddir!=pd.NA)
    
    db_flats = nb_values(db1, max_cols=max(db["col"]), col = ["seqno","ddir"], db_sub = db_flats.reset_index(drop=True))
    
    db_flats = db_flats.loc[np.logical_and(np.logical_and(~pd.isna(db_flats["ddir_n"]),db_flats["ddir_n"]==5),db_flats["n"]!=5)].reset_index(drop=True)
    db_flats["patch"] = pd.NA
    db_flats = db_flats.sort_values(by=["seqno"]).reset_index(drop=True)
    
    

    
    
    if not db_flats.empty:
#         print(db_flats)
        p_n = 1
        for i in range(len(db_flats)):
#             print(db_flats.iloc[i])
            cell = db_flats.iloc[i]["seqno"]
            
            if pd.isna(db_flats.iloc[i]["patch"]):
                db_flats.loc[db_flats["seqno"]==cell,"patch"] = p_n
                p_n += 1
            
            p = db_flats.iloc[i]["patch"] 
            cells_n = np.unique(db_flats.loc[datar.base.testing.is_in(db_flats["seqno_n"],cell),"seqno"])
            p_old = np.unique(db_flats.loc[datar.base.testing.is_in(db_flats["seqno"],cells_n),"patch"])
            p_old = [x for x in p_old if (not x is pd.NA) and (x!=p)]
#             print(p_old)
            db_flats.loc[datar.base.testing.is_in(db_flats["seqno"],cells_n),"patch"] = p
            if len(p_old)>0:
                db_flats.loc[datar.base.testing.is_in(db_flats["patch"],p_old),"patch"] = p
            
            
        # Get middle cell in a patch
        pit_centres = db_flats >> dplyr.select(~f.n,~f.ddir_n,~f.seqno_n)
        pit_centres = pit_centres >> dplyr.distinct()
        

        #POSSIBLE PARALLEL STUFF??? at calculating the root mean square for the distance of each row (point)
        pit_centres = pit_centres >> dplyr.group_by(f.patch)
        
        pit_centres = pit_centres >> dplyr.mutate(centre = [tuple([datar.base.round(datar.base.median(f.col)),datar.base.round(datar.base.median(f.row))])])
        pit_centres = pit_centres >> dplyr.ungroup()
        
        pit_centres["cell"] = pit_centres.apply(lambda x: tuple([x["col"],x["row"]]),axis=1)
        pit_centres["dist"] = pit_centres.apply(lambda x: np.sqrt((x["centre"][0]-x["cell"][0])**2+(x["centre"][1]-x["cell"][1])**2),axis=1)
        pit_centres["dist_min"] = np.amin(pit_centres["dist"])
        pit_centres["n_p"] = len(pit_centres["seqno"])
        
        pit_centres = pit_centres >> dplyr.group_by(f.patch)
        pit_centres = pit_centres >> dplyr.summarise(seqno = f.seqno[f.dist==f.dist_min], n_p = f.n_p.unique())
        
        pit_centres = pit_centres >> dplyr.filter(f.n_p>1)
        
        #MAYBE THIS WILL RECQUIRE A DIFFERENT HANDLING IF THE seqno HAS MORE THAN ONE ELEMENTS. 
        pit_centres = pit_centres["seqno"].tolist()
        
        # Recalculate flow directions for all flat cells in a patch towards the center
        # Iterate over all new flow flat cells until finished
        
        
        # Get directions to pit centers by pit patch
        db_flats = db_flats >> dplyr.group_by(f.seqno)
        db_flats = db_flats >> dplyr.mutate(ddir_opts = [tuple([f.n])])
        db_flats["ddir_opts"] = db_flats["ddir_opts"].apply(lambda x: [y for y in x[0]])
        db_flats = db_flats >> dplyr.select(~f.n,~f.ddir_n,~f.seqno_n)
        db_flats = db_flats >> dplyr.distinct(f.seqno,_keep_all=True)
        db_flats = db_flats >> dplyr.mutate(centre = datar.base.testing.is_in(f.seqno,pit_centres))
        db_flats = db_flats >> dplyr.group_by(f.patch)
        db_flats = db_flats >> dplyr.mutate(row_f=f.row[f.centre].iloc[0], col_f=f.col[f.centre].iloc[0])
        db_flats = db_flats >> dplyr.ungroup()
        
        db_flats["ddir_new"] = db_flats[["row","col","row_f","col_f","ddir_opts"]].apply(lambda x: get_dir(x["row"],x["col"],x["row_f"],x["col_f"],x["ddir_opts"]),axis=1)
        
        db1.loc[(db_flats["seqno"]-1).tolist(),"ddir"] = db_flats["ddir_new"].values
        
        
        

    # Get flow direction (seqno of next cell)
    db1 = flow_values(db1,max_cols=max(db["col"]), col = ["seqno"])
    
    
    
    db1 = db1.rename(columns={"seqno_next": "drec"})
    
    # Fix circular flow among flat cells
    # Shouldn't be necessary anymore, as specified lowest cell already?

    # Check for side-by-side pits, replace so one flows into the other
    # Shouldn't be necessary anymore, as specified lowest cell already?
    
    return db1


