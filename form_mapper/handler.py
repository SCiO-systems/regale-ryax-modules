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
from LITAP_read_write import *
from form_01_calc_form import *
from form_02_calc_weti import *
from form_03_calc_relz import *
from form_04_calc_length import *



def handle(module_input):

    arg1 = "grid"
    arg2 = "str_val"
    arg3 = "ridge_val" 
    arg4 = "output folder"

    output_backup_folder = "../python_outputs/backup/"
    output_stats_folder = "../python_outputs/flow/"

    file = module_input["file"]     #"../../landmapr/LITAP/inst/extdata/testELEV.dbf"
    grid = module_input["grid"]
    str_val = module_input["str_val"]
    ridge_val = module_input["ridge_val"]
    verbose = module_input["verbose"]
    resume = None


    if (resume==None): resume=""

    # # in order to be compatible with ryax platform, inputs with value -1 are translated to None value
    # if clim==-1:
    #     clim = None
    # if rlim==-1:
    #     rlim = None


    #--------------------------------------------------------------------------------------------

    #Load File
    # loading data from flow execution
    # Get backup fill dem
    db = get_previous(folder,step="fill",where="flow")
    db = dplyr.select(db, f.seqno, f.row, f.col, f.elev, f.drec, f.upslope, f.fill_shed, f.local_shed)
    db = add_buffer(db)

    #--------------------------------------------------------------------------------------------

    # Get backup inverted dem
    idb = get_previous(folder,step="ilocal",where="flow")
    if "ldir" in idb.columns:
        idb = datar.all.rename(idb, ddir="ldir")
    idb = dplyr.select(idb,f.seqno,f.row,f.col,f.elev,f.drec,f.ddir,f.upslope,f.shedno)
    idb = add_buffer(idb)

    #--------------------------------------------------------------------------------------------

    # Get backup pond stats
    pond = pd.read_csv(folder + "flow/stats_pond.csv")
    pond = add_buffer(db,stats=pond)


    os.system("mkdir " + folder + "form/")
    #%% # Form 
    if (resume=="" or resume=="form"):
        db_form = calc_form(db, grid,verbose=verbose)
        
        save_output2(data=db_form, name="form", locs=folder, out_format=out_format, where = "form")


    #--------------------------------------------------------------------------------------------

    #%% # Wetness indices 
    if (resume=="" or resume=="weti"):
        db_weti = calc_weti(db, grid, verbose = verbose)

        db_form = dplyr.full_join(db_form, db_weti, by=["seqno", "col", "row", "buffer"])
        
        db_form["lnqarea1"] = np.where(db_form["aspect"] > -1, np.log(db_form["qarea1"].astype(float)), 0)
        db_form["lnqarea2"] = np.where(db_form["aspect"] > -1, np.log(db_form["qarea2"].astype(float)), 0)
        db_form["new_asp"] = np.where(db_form["aspect"] > -1, db_form["aspect"] + 45, 0)
        db_form["new_asp"] = np.where(db_form["new_asp"] > 360,db_form["new_asp"] -360, db_form["new_asp"])
        db_form["lnqarea1"] = round(db_form["lnqarea1"], 3)
        db_form["lnqarea2"] = round(db_form["lnqarea2"], 3)

        #CHECK SAVED FILE THAT IS SIMILAR TO THE R-PRODUCED FILE
        save_output2(data=db_form, name="weti", locs=folder, out_format=out_format, where = "form")
        
        del db_form
        del db_weti

    #--------------------------------------------------------------------------------------------

    #%% # Relief
    if (resume=="" or resume=="relief"):
        db_relz = calc_relz(db, idb, str_val = str_val, ridge_val = ridge_val, pond = pond, verbose = verbose)
        
        save_output2(data=db_relz, name="relief", locs=folder, out_format=out_format, where = "form")
            


    #--------------------------------------------------------------------------------------------

    #%% # Length 
    if (resume=="" or resume=="length"):
        db_length = calc_length(db, db_relz, verbose = verbose)

        save_output2(data=db_length, name="length", locs=folder, out_format=out_format, where = "form")

        del db_length
        del db_relz
         
    successful = 200
    
    return {'successful' : successful}
    