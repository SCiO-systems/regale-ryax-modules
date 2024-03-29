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
from datar import dplyr,tidyr
from datar.all import f
from itertools import compress
import math
import datar


def check_rules(arule,crule):
    arule_weti_relief = ["prof", "plan", "slope_pct", "aspect", "qweti", "qarea", "lnqarea", "new_asp", "pctz2st", "pctz2pit", "z2pit"]
    
    checkif = all([item in arule_weti_relief for item in arule["attr_in"]])
    if not checkif:
        invalid_attr = [not item in arule_weti_relief for item in arule["attr_in"]]
        print("Invalid 'arule' attribute(s): " + arule["attr_in"].iloc[invalid_attr].values)
        raise("Invalid 'arule' attribute(s): " )
        
        
    checkif = all([item in arule["class_out"].values for item in crule["fuzattr"].values])
    if not checkif:
        invalid_attr = [item in arule["class_out"].values for item in crule["fuzattr"].values]
        print("Some arule classes ('class_out') don't match crule fuzzy attributes: " + crule["fuzattr"].iloc[invalid_attr].values)
        print("Arule classes available: ", arule["class_out"].values)
        raise("Invalid 'arule' attribute(s): " )
    
    return