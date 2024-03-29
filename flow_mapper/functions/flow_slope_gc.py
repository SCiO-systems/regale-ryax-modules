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



#flow_slope_gc.R functions
def slope_gc(db, grid=1):
    
    if not datar.base.testing.is_in("buffer",db.columns):
        db = add_buffer(db)
    
    # Calculate Gradients and Curvatures
    db_slopes = nb_values(db, max_cols=max(db["col"]), format_="wide" )
    
    db_slopes.loc[pd.isna(db_slopes["elev_n2"]),"elev_n2"] = db_slopes["elev_n5"]
    db_slopes.loc[pd.isna(db_slopes["elev_n4"]),"elev_n4"] = db_slopes["elev_n5"]
    db_slopes.loc[pd.isna(db_slopes["elev_n6"]),"elev_n6"] = db_slopes["elev_n5"]
    db_slopes.loc[pd.isna(db_slopes["elev_n8"]),"elev_n8"] = db_slopes["elev_n5"]
    
    db_slopes["sgre"] = (db_slopes["elev_n4"]-db_slopes["elev_n6"])/(2*grid)  # Gradient towards east
    db_slopes["sgr"] = abs(db_slopes["sgre"])
    
    db_slopes["sgcn"] = (db_slopes["elev_n2"]-db_slopes["elev_n8"])/(2*grid)  # Gradient towards north
    db_slopes["sgc"] = abs(db_slopes["sgcn"])
    
    db_slopes["scr"] = ((2*db_slopes["elev"] - db_slopes["elev_n4"] - db_slopes["elev_n6"])/grid**2).replace(np.nan,pd.NA)
    db_slopes["scc"] = ((2*db_slopes["elev"] - db_slopes["elev_n2"] - db_slopes["elev_n8"])/grid**2).replace(np.nan,pd.NA)
    
    db_slopes = dplyr.select(db_slopes, ~dplyr.contains("_n"))
    
    # Fix zero values for sgr and sgc and start for sgre and sgcn
    db_slopes.loc[db_slopes["sgr"]==0,"sgr"] = 0.00001
    db_slopes.loc[db_slopes["sgc"]==0,"sgc"] = 0.00001
    
    # Fix sign on small values for sgre - iterate until all solved
    

    n_sgre = np.sum(db_slopes["sgre"]==0)    
    
    while(n_sgre>0):
        db_slopes = nb_values(db_slopes, max_cols=max(db["col"]),col=["sgre"], format_="wide" )
        # If n4 0, iterate over
        db_slopes = db_slopes.replace(np.nan,pd.NA)
        db_slopes["s"] = db_slopes.apply(lambda x: np.where(pd.isna(x["sgre_n4"]),np.sign(x["sgre_n6"]),np.sign(x["sgre_n4"])),axis=1)
        db_slopes.loc[db_slopes["sgre"]==0,"sgre"] = 0.00001 * db_slopes['s']
        db_slopes = db_slopes >> dplyr.select(~dplyr.contains("_n"),~f.s)
        my_sum = np.sum(db_slopes["sgre"]==0)
        if n_sgre == my_sum:
            db_slopes.loc[db_slopes["sgre"]==0,"sgre"] = 0.00001
            n_sgre = 0
        else:
            n_sgre = my_sum
    
    # Fix zero values for sgcn - iterate until all solved
    n_sgcn = np.sum(db_slopes["sgcn"]==0)    
    
    while(n_sgcn>0):
        db_slopes = nb_values(db_slopes, max_cols=max(db["col"]),col=["sgcn"], format_="wide" )
        # If n4 0, iterate over
        db_slopes = db_slopes.replace(np.nan,pd.NA)
        db_slopes["s"] = db_slopes.apply(lambda x: np.where(pd.isna(x["sgcn_n2"]),np.sign(x["sgcn_n8"]),np.sign(x["sgcn_n2"])),axis=1)
        db_slopes.loc[db_slopes["sgcn"]==0,"sgcn"] = 0.00001 * db_slopes['s']
        db_slopes = db_slopes >> dplyr.select(~dplyr.contains("_n"),~f.s)
        my_sum = np.sum(db_slopes["sgcn"]==0)
        if n_sgcn == my_sum:
            db_slopes.loc[db_slopes["sgcn"]==0,"sgcn"] = 0.00001
            n_sgcn = 0
        else:
            n_sgcn = my_sum
        
        
    # Hillslope records
    # Get east/west or north/south facing
    db_slopes = db_slopes >> dplyr.mutate(hill_r_dir = dplyr.if_else(f.sgre>0,2,4),  # 2 = east, 4 = west
                                          hill_c_dir = dplyr.if_else(f.sgcn>0,1,3))  # 1 = north, 3 = south
    
    db_slopes.loc[pd.isna(db_slopes["elev"]),"hill_r_n"] = pd.NA    
    db_slopes.loc[pd.isna(db_slopes["elev"]),"hill_r_n"] = pd.NA
    
    db_slopes["sgre_sign"] = db_slopes["sgre"].apply(lambda x: np.sign(x))
    db_slopes["sgcn_sign"] = db_slopes["sgcn"].apply(lambda x: np.sign(x))
    
    
    # Label east/west hillslopes
    db_slopes = db_slopes >> dplyr.group_by(f.row)
    
    #substracting 1 so as to achieve the same results as the R code, 
    #TO CHECK if something is wrong later on
    db_slopes = db_slopes >> dplyr.mutate(lag_s = dplyr.lag(f.sgre_sign),
                                         hill_r_n = np.cumsum(f.sgre_sign!=f.lag_s or pd.isna(f.lag_s))-1)
    
    db_slopes.loc[pd.isna(db_slopes["elev"]),"hill_r_n"] = pd.NA
    
    # Label east/west cells within a hillslope
    
    db_slopes = db_slopes >> dplyr.ungroup()

    db_slopes["hill_r_cell"] = db_slopes.groupby(["row","hill_r_n"]).cumcount() + 1
    

    
    
    # Label north/south hillslopes
    db_slopes = db_slopes >> dplyr.group_by(f.col)
    
    #substracting 1 so as to achieve the same results as the R code, 
    #TO CHECK if something is wrong later on
    db_slopes = db_slopes >> dplyr.mutate(lag_s = dplyr.lag(f.sgcn_sign),
                                         hill_c_n = np.cumsum(f.sgcn_sign!=f.lag_s or pd.isna(f.lag_s))-1)
    
    db_slopes.loc[pd.isna(db_slopes["elev"]),"hill_c_n"] = pd.NA
    
    # Label east/west cells within a hillslope
    
    db_slopes = db_slopes >> dplyr.ungroup()

    db_slopes["hill_c_cell"] = db_slopes.groupby(["col","hill_c_n"]).cumcount() + 1
    
    # Remove labels on missing values
    db_slopes.loc[pd.isna(db_slopes["elev"]),"hill_r_cell"] = pd.NA
    db_slopes.loc[pd.isna(db_slopes["elev"]),"hill_c_cell"] = pd.NA
    
    db_slopes = db_slopes >> dplyr.select(~f["lag_s","sgre_sign","sgcn_sign"])
    
    return db_slopes