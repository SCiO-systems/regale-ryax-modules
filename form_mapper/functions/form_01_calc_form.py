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


def rad_deg(x):
        return x*180/math.pi

def deg_rad(x):
    return x*math.pi/180

def calc_form(db,grid=10,verbose=False):
    
    # slope, aspect, prof, plan 
    db = dplyr.select(db, f.seqno,f.elev,f.row,f.col,f.buffer)
    db = nb_values(db, max_cols=np.max(db["col"]),format_="wide")
    
    for column in db.columns:
        if "elev_n" in column:
            db[column] = 100 * db[column]
    #next R command calculates in parallel how many nan values are in elev_n* columns
    elev_columns = [col for col in db.columns if 'elev_n' in col]
    sum_elev = db[elev_columns].replace(np.nan,pd.NA)
    sum_elev = sum_elev.isnull()
    sum_elev = sum_elev.sum(axis=1)
    db["sum_elev"] = sum_elev
    
    db["slope_x"] = (db["elev_n6"] - db["elev_n4"])/(2*grid)
    db["slope_y"] = (db["elev_n2"] - db["elev_n8"])/(2*grid)
    db["slope_pct"] = np.sqrt(db["slope_x"]**2 + db["slope_y"]**2)
    db["slope_deg"] = rad_deg(np.arctan(db["slope_pct"]/100))
    db["aspect"] = aspect(db["slope_x"],db["slope_y"],db["slope_pct"])
    
    db["prof_aspect"] = dplyr.if_else(db["aspect"]> 180, db["aspect"] - 180, db["aspect"])
    
    db["plan_aspect"] = dplyr.if_else((db["prof_aspect"] + 90) > 180, db["prof_aspect"] + 90 - 180, db["prof_aspect"] + 90)
    
    db["prof"] = dplyr.if_else(db["sum_elev"] > 0, pd.NA, prof_plan(db["prof_aspect"], db["elev_n1"], db["elev_n2"],
                                                  db["elev_n3"], db["elev_n4"], db["elev_n5"],
                                                  db["elev_n6"], db["elev_n7"], db["elev_n8"],
                                                  db["elev_n9"], grid, db["slope_pct"]))
    db["plan"] = dplyr.if_else(db["sum_elev"] > 0, pd.NA, prof_plan(db["plan_aspect"], db["elev_n1"], db["elev_n2"],
                                                  db["elev_n3"], db["elev_n4"], db["elev_n5"],
                                                  db["elev_n6"], db["elev_n7"], db["elev_n8"],
                                                  db["elev_n9"], grid, db["slope_pct"]))

    db = dplyr.select(db, f.seqno,f.row,f.col,f.slope_pct,f.slope_deg,f.aspect,f.prof,f.plan,f.buffer)
    
    db = db.replace(pd.NA,np.nan)
    
    db["slope_pct"] = db["slope_pct"].round(3)
    db["slope_deg"] = db["slope_deg"].round(3)
    db["aspect"] = db["aspect"].apply(np.round)
    db["prof"] = db["prof"].round(3)
    db["plan"] = db["plan"].round(3)
    
    
    db = db.replace(np.nan,pd.NA)
    
    # First/last rows and cols get adjacent values
    vals = ["slope_pct", "slope_deg", "aspect", "prof", "plan"]
    
    # Note that first and last row over write corners (assume buffer)
    
    # Left Column
    db.loc[db["col"]==2,vals] = db.loc[db["col"]==3,vals].set_index(db.loc[db["col"]==2,vals].index) 
    # Right Column
    db.loc[db["col"]==max(db["col"])-1,vals] = db.loc[db["col"]==max(db["col"])-2,vals].set_index(db.loc[db["col"]==max(db["col"])-1,vals].index)
    # First row
    db.loc[db["row"]==2,vals] = db.loc[db["row"]==3,vals].set_index(db.loc[db["row"]==2,vals].index)
    # Last row
    db.loc[db["row"]==max(db["row"])-1,vals] = db.loc[db["row"]==max(db["row"])-2,vals].set_index(db.loc[db["row"]==max(db["row"])-1,vals].index)
    
    return db

def prof_plan(aspect, elev_n1, elev_n2, elev_n3,
              elev_n4, elev_n5, elev_n6, elev_n7,
              elev_n8, elev_n9, grid, slope_pct):
    if aspect.dtype=="object":
        aspect = pd.to_numeric(aspect)
    x1 = 2 + np.sin(deg_rad(aspect))
    y1 = 2 - np.cos(deg_rad(aspect))
    x2 = 2 - np.sin(deg_rad(aspect))
    y2 = 2 + np.cos(deg_rad(aspect))
    
    z1 = dplyr.case_when(
        aspect<=90, ((2 - y1) * ((elev_n9 * (x1 - 2)) + (elev_n8 * (3 - x1)))) + 
        ((y1 - 1) * ((elev_n6 * (x1 - 2)) + (elev_n5 * (3 - x1)))),
        aspect>90, ((3 - y1) * ((elev_n6 * (x1 - 2)) + (elev_n5 * (3 - x1)))) + 
        ((y1 - 2) * ((elev_n3 * (x1 - 2)) + (elev_n2 * (3 - x1))))
    )
    
    z2 = dplyr.case_when(  
        aspect <= 90, ((3 - y2) * ((elev_n5 * (x2 - 1)) + (elev_n4 * (2 - x2)))) +
        ((y2 - 2) * ((elev_n2 * (x2 - 1)) + (elev_n1 * (2 - x2)))),
        aspect > 90, ((2 - y2) * ((elev_n8 * (x2 - 1)) + (elev_n7 * (2 - x2)))) +
        ((y2 - 1) * ((elev_n5 * (x2 - 1)) + (elev_n4 * (2 - x2))))
    )
    
    
    return dplyr.case_when(slope_pct <= 0, 0,
                           True, rad_deg(np.arctan(((2 * elev_n5) - z1 - z2) / (grid * grid))))

def aspect(slope_x, slope_y, slope_pct):
    local_angle = rad_deg(np.arccos(abs(slope_x)/slope_pct))
    return dplyr.case_when(slope_pct <= 0, 360,
                           np.logical_and(slope_x > 0,slope_y > 0), 270 + local_angle,
                           np.logical_and(slope_x > 0, slope_y < 0), 270 - local_angle,
                           np.logical_and(slope_x < 0, slope_y > 0), 90 - local_angle,
                           np.logical_and(slope_x < 0, slope_y < 0), 90 + local_angle,
                           np.logical_and(slope_x < 0, slope_y == 0), 90,
                           np.logical_and(slope_x > 0, slope_y == 0), 270,
                           np.logical_and(slope_x == 0, slope_y < 0), 180,
                           np.logical_and(slope_x == 0, slope_y > 0), 360,
                           True, pd.NA)