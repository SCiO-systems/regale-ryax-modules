#flow_mapper.py
import sys
import os

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(1, 'functions/')

#imports
from dbfread import DBF
import pandas as pd
import numpy as np
import copy
import string
import datar
from datar import dplyr,tidyr
from datar.all import f
from itertools import compress
import math
import os
import json
import time

#import from other files
from LITAP_functions import *
from LITAP_utils import *
from LITAP_load import *
from LITAP_read_write import *

from flow_calc_ddir import *
from flow_calc_shed import *
from flow_pit_stat_fast import *
from flow_remove_pit_fast import *
from flow_slope_gc import *

import warnings
warnings.filterwarnings("ignore")


#%%
def handle(module_input):
    with open(module_input["input_json"], 'r') as file:
      data = json.load(file)

    file = module_input["input_file"]

    nrow = data["hyperparameters"]["nrow"]
    ncol = data["hyperparameters"]["ncol"]
    missing_value = data["hyperparameters"]["nodata"]
    clim = data["hyperparameters"]["clim"]
    rlim = data["hyperparameters"]["rlim"]
    verbose = data["hyperparameters"]["verbose"]
    verbose = True
    resume = None
    max_area = data["hyperparameters"]["max_area"]
    max_depth = data["hyperparameters"]["max_depth"]

    # TMP_DIR = "/tmp/"
    TMP_DIR = "../data/"

    out_directory = TMP_DIR + "test_dem_outputs/"

    output_flow_folder = out_directory + "flow/"
    output_backup_folder = out_directory + "backup/"

    os.makedirs(output_flow_folder,exist_ok=True)
    os.makedirs(output_backup_folder,exist_ok=True)
    os.makedirs(out_directory + "form/",exist_ok = True)
    os.makedirs(out_directory + "facet/",exist_ok = True)

    if (resume==None): resume=""

    # in order to be compatible with ryax platform, inputs with value -1 are translated to None value
    if clim==-1:
        clim = None
    if rlim==-1:
        rlim = None

    #--------------------------------------------------------------------------------------------        
    ##%% ### Loading File and PreProcessing it

    #Load File
    db_start = load_file(file,nrow=nrow,ncol=ncol,missing_value=missing_value,clim=clim,rlim=rlim,verbose=verbose)
    ncol_orig = ncol
    nrow_orig = nrow
    ncol = max(db_start["col"]) - 2
    nrow = max(db_start["row"]) - 2

    #--------------------------------------------------------------------------------------------        
    #%% Calculate Directions
    print("Calculate Directions")
    if (resume=="" or resume=="directions"):
        db_dir = calc_ddir2(copy.deepcopy(db_start),verbose=verbose)
        
        db_dir.to_csv(output_backup_folder+"dir.csv",index=False)

    #--------------------------------------------------------------------------------------------        
    #%% Calculate watersheds
    print("Calculate watersheds")
    if(resume == "" or resume == "watersheds"):
        db_initial = dict()
        
        db_initial["db"] = calc_shed4(copy.deepcopy(db_dir))

        pit_stats = pit_stat1(copy.deepcopy(db_initial["db"]))
        
        db_initial["stats"] = out_stat(copy.deepcopy(pit_stats))
        # Calc stats for first vol2fl
        db_initial["db"] = calc_vol2fl(copy.deepcopy(db_initial["db"]), i_stats = db_initial["stats"], verbose = verbose)
        
        db_initial["db"].to_csv(output_backup_folder+"initial.csv",index=False)
        db_initial["stats"].to_csv(output_flow_folder+"stats_initial.csv",index=False)

    #--------------------------------------------------------------------------------------------        
    #%% Remove initial pits
    print("Remove initial pits")
    if(resume == "" or resume == "local"):
    #     print(db_initial["db"].iloc[4400:4404])
        db_local = first_pitr1(copy.deepcopy(db_initial["db"]), max_area=max_area,max_depth=max_depth,verbose=verbose)
        stats_local = pit_stat1(copy.deepcopy(db_local))
        stats_local = out_stat(copy.deepcopy(stats_local))
        
        db_local = {
            "db" : db_local,
            "stats" : stats_local
        }
        
        db_local["db"].to_csv(output_backup_folder+"local.csv",index=False)
        db_local["stats"].to_csv(output_backup_folder+"stats_local.csv",index=False)
        
    #--------------------------------------------------------------------------------------------        
    #%% Calc pond Sheds
    print("Calc pond Sheds")
    if(resume == "" or resume == "pond"):
        
        if (len(np.unique(db_local["db"]["shedno"][~pd.isna(db_local["db"]["shedno"])])) >1):
            db_pond = second_pitr1(copy.deepcopy(db_local["db"]),verbose=verbose)
        else:
            db_pond = dict()
            db_pond["db"] = db_local["db"] >> dplyr.mutate(pond_shed = f.local_shed)
            db_pond["stats"] = None
            raise("TODO: check if this is correct")
            
            
        db_pond["db"].to_csv(output_backup_folder+"pond.csv",index=False)
        db_pond["stats"].to_csv(output_backup_folder+"stats_pond.csv",index=False)
        
    #--------------------------------------------------------------------------------------------        
    #%% Calc fill Sheds
    print("Calc fill Sheds")
    if(resume == "" or resume == "fill"):
        
        if (len(np.unique(db_local["db"]["shedno"][~pd.isna(db_local["db"]["shedno"])])) >1):
            # Add pond sheds details to local sheds
            db_local["db"][["vol2fl", "mm2fl", "parea"]] = db_pond["db"][["vol2fl", "mm2fl", "parea"]]
            db_local["db"]["pond_shed"] = db_pond["db"]["pond_shed"]
            
            db_fill = third_pitr1(copy.deepcopy(db_local["db"]),verbose=verbose)
        else:
            print("  Only a single watershed: No fill outputs")
            db_fill = dict()
            db_fill["db"] = db_pond["db"] >> dplyr.mutate(fill_shed = f.local_shed, vol2fl = 0, mm2fl = 0, parea = 0)
            db_fill["stats"] = pd.DataFrame()
         
            
        # Calculate slope gradients and curvatures
        db_fill["db"] = slope_gc(copy.deepcopy(db_fill["db"]), grid=1)
        
        db_fill["db"].to_csv(output_backup_folder+"fill.csv",index=False)
        
        db_fill["stats"].to_csv(output_backup_folder+"stats_fill.csv",index=False)

        #saving pit file TODO the savings and readings at some point
        if len(db_fill["stats"])>0:
            # Create PIT file
            pit = db_fill["stats"]
            pit = datar.all.filter(pit,f.final==True)
            pit["edge_pit"] = False
            pit = datar.all.arrange(pit,f.shedno)
              
            db_fill["db"].to_csv(output_backup_folder+"pit.csv",index=False)
            # pit.to_csv(output_flow_folder+"stats_pit.csv",index=False)
            pit.to_csv(output_backup_folder+"stats_pit.csv",index=False)

    #--------------------------------------------------------------------------------------------        
    #%% ### Inverted DEM
    print("Inverted DEM")
    def invert(db):
        max_elev = np.amax(db["elev"])
        db["elev"] = max_elev - db["elev"]
        return db

    if(resume == "" or resume == "inverted"):
        db_invert = db_local["db"][["elev", "seqno", "row", "col", "missing", "buffer", "elev_orig", "edge_map"]]
        db_invert = invert(copy.deepcopy(db_invert))
        
        # Inverted Directions
        db_idir = calc_ddir2(db_invert, verbose=verbose)

        db_idir.to_csv(output_backup_folder+"idir.csv",index=False)
        
    #--------------------------------------------------------------------------------------------                
    #%% Inverted Watersheds
    print("Inverted Watersheds")
    if(resume == "" or resume == "iwatersheds"):
        db_iinitial = calc_shed4(copy.deepcopy(db_idir))
        db_iinitial.to_csv(output_backup_folder+"iinitial.csv",index=False)

    #--------------------------------------------------------------------------------------------        
    #%% Invert Remove Initial Pits
    print("Invert Remove Initial Pits")
    db_ilocal = first_pitr1(copy.deepcopy(db_iinitial),max_area = max_area, max_depth = max_depth, verbose = verbose)
    if (len(np.unique(db_ilocal["shedno"][~pd.isna(db_ilocal["shedno"])])) >1):
        ipit = pit_stat1(db_ilocal)
        ipit = out_stat(ipit)
        ipit = dplyr.mutate(ipit,edge_pit=False)
    else:
        ipit = pd.DataFrame()

  
    db_ilocal.to_csv(output_backup_folder+"ilocal.csv",index=False)
    ipit.to_csv(output_backup_folder+"stats_ilocal.csv",index=False)

    #--------------------------------------------------------------------------------------------        
    #%% Save output in the right format
    save_output(output_backup_folder)
    print("FlowMapR finished execution.")

    return {'python_outputs' : out_directory}


#%%

f1 = {"input_json" :"../data/flow_test_dem_input_json.json",
  "input_file" : "../data/test_dem.tif"
  }
  

start_time = time.time()
t = handle(f1)  
end_time = time.time()
print("Total time:", str(end_time-start_time))