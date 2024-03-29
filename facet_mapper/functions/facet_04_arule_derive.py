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

def arule_derive(weti, relief, n_remove):
    perc = dplyr.left_join(weti, dplyr.select(relief, f.seqno, f.pctz2st, f.pctz2pit, f.z2pit), by="seqno")
    perc = datar.all.filter(perc,f.row>(n_remove+1),f.row<(max(perc["row"]) - n_remove-1), f.col>(n_remove+1),f.col<(max(perc["col"]) - n_remove-1)).reset_index(drop=True)
    
    perc = dplyr.select(perc, f.prof, f.plan, f.slope_pct, f.qweti, f.pctz2st, f.pctz2pit, f.z2pit)
    perc = perc.rename(columns={"slope_pct" : "slope"})
    # dplyr.summarize code // from here
    perc["n"] = len(perc)
    t = datar.base.arithmetic.quantile(perc,probs=[0.1,0.25,0.35,0.50,0.65,0.70,0.75,0.90]).reset_index(drop=True)    
    quantile_list = ["p10","p25","p35","p50","p65","p70","p75","p90"]
    df = {"n" : len(perc)}
    for col in t.columns:
        
        for i,quantile in enumerate(quantile_list):
            temp = col + "_" + quantile
            df[temp] = [t[col].iloc[i]]
    
    perc = pd.DataFrame.from_dict(df)
    # to here //
    
    # arule_template code // from here
    arule_template = {
        "sortorder" : list(range(1,18)),
        "file_in" : ["formfile","formfile","formfile","formfile","formfile","formfile","formfile","formfile","formfile","formfile","relzfile","relzfile","relzfile","relzfile","relzfile","relzfile","relzfile"],
        "attr_in" : ["PROF","PROF","PROF","PLAN","PLAN","PLAN","QWETI","QWETI","SLOPE","SLOPE","PCTZ2ST","PCTZ2ST","PCTZ2ST","PCTZ2PIT","PCTZ2PIT","PCTZ2PIT","Z2PIT"],
        "class_out" : ["CONVEX_D","CONCAVE_D","PLANAR_D","CONVEX_A","CONCAVE_A","PLANAR_A","HIGH_WI","LOW_WI","NEAR_LEVEL","REL_STEEP","NEAR_DIV","NEAR_HALF","NEAR_CHAN","NEAR_PEAK","NEAR_MID","NEAR_PIT","HI_ABOVE"],
        "model" : [4,5,1,4,5,1,4,5,5,4,4,1,5,4,1,5,4],
        "calc" : ["bd1","bd2","lhd","bd1","bd2","lhd","bd1","bd2","bd2","bd1","bd1","lhd","bd2","bd1","lhd","bd2","bd1"]
        }
    a = pd.DataFrame.from_dict(arule_template)
    # to here //
    
    def big_or_min(val,cuttoff):
        return dplyr.if_else(abs(val)<abs(cuttoff), cuttoff, val)[0]
            
    b = [
        big_or_min(perc["prof_p90"], 0.1),
        big_or_min(perc["prof_p10"],-0.1),
        0,
        big_or_min(perc["plan_p90"], 0.1),
        big_or_min(perc["plan_p10"],-0.1),
        0,
        perc["qweti_p90"].values[0],
        perc["qweti_p10"].values[0],
        0,
        perc["slope_p90"].values[0],
        100,
        50,
        0,
        100,
        50,
        0,
        perc["z2pit_p90"].values[0]
        ]
    b_low = [0,0,0,0,0,0,0,0,0,0,0,50,0,0,50,0,0]
    b_hi = [0,0,0,0,0,0,0,0,0,0,0,50,0,0,50,0,0]
    
    d = [
        big_or_min((perc["prof_p90"].values[0] - perc["prof_p65"].values[0])/2, 0.01),
        big_or_min((perc["prof_p35"].values[0] - perc["prof_p10"].values[0])/2, 0.01),
        big_or_min((perc["prof_p75"].values[0] - perc["prof_p25"].values[0])/2, 0.01),
        big_or_min((perc["plan_p90"].values[0] - perc["plan_p65"].values[0])/2, 0.01),
        big_or_min((perc["plan_p35"].values[0] - perc["plan_p10"].values[0])/2, 0.01),
        big_or_min((perc["plan_p75"].values[0] - perc["plan_p25"].values[0])/2, 0.01),
        (perc["qweti_p90"].values[0] - perc["qweti_p50"].values[0])/2,
        (perc["qweti_p50"].values[0] - perc["qweti_p10"].values[0])/2,
        big_or_min(perc["slope_p25"].values[0]/2, 0.01),
        big_or_min((perc["slope_p90"].values[0] - perc["slope_p25"].values[0])/2, 0.01),
        (100 - perc["pctz2st_p75"].values[0]) / 2,
        (perc["pctz2st_p75"].values[0] - perc["pctz2st_p25"].values[0]) / 2,
        perc["pctz2st_p25"].values[0] / 2,
        (100 - perc["pctz2pit_p75"].values[0]) / 2,
        (perc["pctz2pit_p75"].values[0] - perc["pctz2pit_p25"].values[0]) / 2,
        perc["pctz2pit_p25"].values[0] / 2,
        (perc["z2pit_p90"].values[0] - perc["z2pit_p70"].values[0]) / 2        
        ]
    

    a["b"] = b
    a["b_low"] = b_low
    a["b_hi"] = b_hi
    a["d"] = d
    
    def b_calcs(calc, b, d, b_low, b_hi, btype):
        if btype == 1:
            x = dplyr.case_when(calc == "bd1", b - d,
                                calc == "bd2", 0,
                                calc == "lhd", b_low - d)
        if btype == 2:
            x = dplyr.case_when(calc == "bd1", 0,
                                calc == "bd2", b + d,
                                calc == "lhd", b_hi + d)
        return x
    
    b1 = b_calcs(a["calc"],a["b"],a["d"], a["b_low"],a["b_hi"],1)
    b2 = b_calcs(a["calc"],a["b"],a["d"], a["b_low"],a["b_hi"],2)
    
    a["b1"] = b1
    a["b2"] = b2
    
    a = a.drop('calc', 1)
    
    return a
