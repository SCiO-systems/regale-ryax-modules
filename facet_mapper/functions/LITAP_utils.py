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
from LITAP_load import *

# LITAP_utils.R functions
def na_omit(x):
    return x.dropna()

    
def max_na(x):
    
    if sum(np.isnan(x))==len(x):
        y = np.nan
    else:
        y = np.nanmax(x)    
    
    return y