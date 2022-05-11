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
from LITAP_read_write import *
from LITAP_functions import *




def calc_length(db,relz, grid=5, verbose=False):
    relz = dplyr.left_join(relz, db[["seqno","row","col"]], by = ["seqno","row","col"])
    
    
    relz["l2pit"] = datar.all.sqrt((relz["pit_col"] - relz["col"])**2 + (relz["pit_row"] - relz["row"])**2) * grid
    relz["l2pit"] = datar.all.sqrt(relz["l2pit"]**2 + (relz["z2pit"] * grid)**2)
    
    relz["l2peak"] = datar.all.sqrt((relz["peak_col"] - relz["col"])**2 + (relz["peak_row"] - relz["row"])**2) * grid
    relz["l2peak"] = datar.all.sqrt(relz["l2peak"]**2 + (relz["z2peak"] * grid)**2)
    
    relz["l2str"] = datar.all.sqrt((relz["str_col"] - relz["col"])**2 + (relz["str_row"] - relz["row"])**2) * grid
    relz["l2str"] = datar.all.sqrt(relz["l2str"]**2 + (relz["z2st"] * grid)**2)
    
    relz["l2div"] = datar.all.sqrt((relz["cr_col"] - relz["col"])**2 + (relz["cr_row"] - relz["row"])**2) * grid
    relz["l2div"] = datar.all.sqrt(relz["l2div"]**2 + (relz["z2cr"] * grid)**2)
    # In C++ version, z2peak, but in foxitpro, z2cr. USE z2cr which I think is more correct
    
    relz["lpit2peak"] = relz["l2pit"] + relz["l2peak"]
    relz["lstr2div"] = relz["l2str"] + relz["l2div"]
    
    relz["ppit2peakl"] = datar.all.trunc((relz["l2pit"] / relz["lpit2peak"]) * 100)
    relz["pstr2divl"] = datar.all.trunc((relz["l2str"] / relz["lstr2div"]) * 100)
    
    
    relz["l2pit"] = relz["l2pit"].round(1)
    relz["l2peak"] = relz["l2peak"].round(1)
    relz["lpit2peak"] = relz["lpit2peak"].round(1)
    relz["l2str"] = relz["l2str"].round(1)
    relz["l2div"] = relz["l2div"].round(1)
    relz["lstr2div"] = relz["lstr2div"].round(1)
    
    relz.loc[relz["lpit2peak"]==0,"ppit2peakl"] = 0
    relz.loc[relz["lstr2div"]==0,"pstr2divl"] = 0

    return relz