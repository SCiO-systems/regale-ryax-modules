#flow_mapper.py

import sys
sys.path.insert(1, './functions/')

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
from flow_calc_ddir import *
from flow_calc_shed import *
from flow_pit_stat_fast import *
from flow_remove_pit_fast import *
from flow_slope_gc import *


def handle(module_input):

    arg1 = "DEM file"
    arg2 = "nrows"
    arg3 = "ncols"
    arg4 = "nodata"
    arg5 = "max_area"
    arg6 = "max_depth"
    arg7 = "output folder"


    output_backup_folder = "../python_outputs/backup/"
    output_stats_folder = "../python_outputs/flow_stats/"

    file = module_input["file"]     #"../../landmapr/LITAP/inst/extdata/testELEV.dbf"
    nrow = module_input["nrow"]
    ncol = module_input["ncol"]
    missing_value = module_input["missing_value"]
    clim = module_input["clim"]
    rlim = module_input["rlim"]
    verbose = module_input["verbose"]
    resume = None
    max_area = module_input["max_area"]
    max_depth = module_input["max_depth"]


    if (resume==None): resume=""

    # in order to be compatible with ryax platform, inputs with value -1 are translated to None value
    if clim==-1:
        clim = None
    if rlim==-1:
        rlim = None


    #--------------------------------------------------------------------------------------------

    #Load File
    db_start = load_file(file,nrow=nrow,ncol=ncol,missing_value=missing_value,clim=clim,rlim=rlim,verbose=verbose)
    ncol_orig = ncol
    nrow_orig = nrow
    ncol = max(db_start["col"]) - 2
    nrow = max(db_start["row"]) - 2

    #--------------------------------------------------------------------------------------------

    #Calculate Directions
    if (resume=="" or resume=="directions"):
        db_dir = calc_ddir2(copy.copy(db_start),verbose=verbose)
        db_dir.to_csv(output_backup_folder+"dir.csv",index=False)

    #--------------------------------------------------------------------------------------------

    # Calculate watersheds

    if(resume == "" or resume == "watersheds"):
        db_initial = dict()
        
        db_initial["db"] = calc_shed4(copy.copy(db_dir))
    #     print(db_initial["db"].iloc[4400:4403])
    # db_initial["db"] 
        pit_stats = pit_stat1(copy.copy(db_initial["db"]))
    #     print(db_initial["db"].iloc[4400:4403])
        
        db_initial["stats"] = out_stat(copy.copy(pit_stats))
    #     print(db_initial["stats"])
        # Calc stats for first vol2fl
        db_initial["db"] = calc_vol2fl(copy.copy(db_initial["db"]), i_stats = db_initial["stats"], verbose = verbose)
        
        db_initial["db"].to_csv(output_backup_folder+"initial.csv",index=False)
        db_initial["stats"].to_csv(output_stats_folder+"stats_initial.csv",index=False)
    #     print(db_initial["db"].iloc[4400:4404])    


    #--------------------------------------------------------------------------------------------

    # Remove initial pits
    if(resume == "" or resume == "local"):
    #     print(db_initial["db"].iloc[4400:4404])
        db_local = first_pitr1(copy.copy(db_initial["db"]), max_area=max_area,max_depth=max_depth,verbose=verbose)
        stats_local = pit_stat1(copy.copy(db_local))
        stats_local = out_stat(copy.copy(stats_local))
        
        db_local = {
            "db" : db_local,
            "stats" : stats_local
        }
        
        db_local["db"].to_csv(output_backup_folder+"local.csv",index=False)
        db_local["stats"].to_csv(output_stats_folder+"stats_local.csv",index=False)
    # db_local




    #--------------------------------------------------------------------------------------------

    # Calc pond Sheds
    if(resume == "" or resume == "pond"):
        
        if (len(np.unique(db_local["db"]["shedno"][~pd.isna(db_local["db"]["shedno"])])) >1):
            db_pond = second_pitr1(copy.copy(db_local["db"]),verbose=verbose)
        else:
            db_pond = dict()
            db_pond["db"] = db_local["db"] >> dplyr.mutate(pond_shed = f.local_shed)
            db_pond["stats"] = None
            raise("TODO: check if this is correct")
            
            
        db_pond["db"].to_csv(output_backup_folder+"pond.csv",index=False)
        db_pond["stats"].to_csv(output_stats_folder+"stats_pond.csv",index=False)
            


    #--------------------------------------------------------------------------------------------

    # Calc fill Sheds
    if(resume == "" or resume == "fill"):
        
        if (len(np.unique(db_local["db"]["shedno"][~pd.isna(db_local["db"]["shedno"])])) >1):
            # Add pond sheds details to local sheds
            db_local["db"][["vol2fl", "mm2fl", "parea"]] = db_pond["db"][["vol2fl", "mm2fl", "parea"]]
            db_local["db"]["pond_shed"] = db_pond["db"]["pond_shed"]
            
            db_fill = third_pitr1(copy.copy(db_local["db"]),verbose=verbose)
        else:
            print("  Only a single watershed: No fill outputs")
            db_fill = dict()
            db_fill["db"] = db_pond["db"] >> dplyr.mutate(fill_shed = f.local_shed, vol2fl = 0, mm2fl = 0, parea = 0)
            db_fill["stats"] = pd.DataFrame()
         
            
        # Calculate slope gradients and curvatures
        db_fill["db"] = slope_gc(copy.copy(db_fill["db"]), grid=1)
        
        db_fill["db"].to_csv(output_backup_folder+"fill.csv",index=False)
        db_fill["stats"].to_csv(output_stats_folder+"stats_fill.csv",index=False)       

        #saving pit file TODO the savings and readings at some point
        if len(db_fill["stats"])>0:
            # Create PIT file
            pit = db_fill["stats"]
            pit = pit >> dplyr.filter(f.final==True)
            pit["edge_pit"] = False
            pit = pit >> dplyr.arrange(f.shedno)
              
            db_fill["db"].to_csv(output_backup_folder+"pit.csv",index=False)
            pit.to_csv(output_stats_folder+"stats_pit.csv",index=False)       



    #--------------------------------------------------------------------------------------------

    def invert(db):
        max_elev = np.amax(db["elev"])
        db["elev"] = max_elev - db["elev"]
        return db



    # Inverted DEM 
    if(resume == "" or resume == "inverted"):
        db_invert = db_local["db"][["elev", "seqno", "row", "col", "missing", "buffer", "elev_orig", "edge_map"]]
        db_invert = invert(copy.copy(db_invert))
        
        # Inverted Directions
        db_idir = calc_ddir2(db_invert, verbose=verbose)
        db_idir.to_csv(output_backup_folder+"idir.csv",index=False)
        
        


    #--------------------------------------------------------------------------------------------

    # Inverted Watersheds
    if(resume == "" or resume == "iwatersheds"):
        db_iinitial = calc_shed4(copy.copy(db_idir))
        db_iinitial.to_csv(output_backup_folder+"iinitial.csv",index=False)



    #--------------------------------------------------------------------------------------------


    # Invert Remove Initial Pits
    db_ilocal = first_pitr1(copy.copy(db_iinitial),max_area = max_area, max_depth = max_depth, verbose = verbose)
    if (len(np.unique(db_ilocal["shedno"][~pd.isna(db_ilocal["shedno"])])) >1):
        ipit = pit_stat1(db_ilocal)
        ipit = out_stat(ipit)
        ipit = ipit >> dplyr.mutate(edge_pit=False)
    else:
        ipit = pd.DataFrame()

        
    db_ilocal.to_csv(output_backup_folder+"ilocal.csv",index=False)
    ipit.to_csv(output_stats_folder+"stats_ilocal.csv",index=False)


    #--------------------------------------------------------------------------------------------

    db_fill["db"].to_csv(output_stats_folder+"dem_fill.csv",index=False)
    db_ilocal.to_csv(output_stats_folder+"dem_ilocal.csv",index=False)

    return db_ilocal