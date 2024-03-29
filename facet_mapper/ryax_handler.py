#facet_mapper.py

import sys
import os

# os.chdir("C:/Users/Christos/Desktop/SCiO_REGALE/REGALE/regale-ryax-modules/facet_mapper_version_2/")
# print(os.listdir("./"))
# print(os.listdir("../../final/python_epirus_3_example_outputs/"))
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(1, 'functions/')

#imports
import pandas as pd
import numpy as np
from datar import dplyr,tidyr
from datar.all import f
import datar
from tifffile import imsave,imwrite
import json
import time

#import from other files
from LITAP_functions import *
from LITAP_utils import *
from LITAP_load import *
from LITAP_read_write import *
from facet_01_lsm import *
from facet_03_check_facet import *
from facet_04_arule_derive import *


import warnings
warnings.filterwarnings("ignore")

#%%
def handle(module_input):

    # Initializing parameters

    # folder Character. Location of [facet_mapper()] output
    # arule Character. Location of ARULE file. If NULL, A Rules are derived from the dem file (see Details).
    # crule Character. Location of CRULE file
    # n_remove Numeric. 
    # Number of cells (rows/columns) to remove around the edge of the dem before deriving the A Rules.
    # procedure Character. Which LSM procedure to use. One of `lsm` (Original LandMapR program) or `bc_pem` (newer BC-PEM Direct-to-Site-SEries DSS program).
    # ARULE zone file. If `procedure = "bc_pem"`, zones must either be defined for each seqno in the weti dem file, OR must be provided as an index file here. 
    # With a `zone` defined for each `seqno`. `zone` file can be either dem (.dem), Excel (.xlsx, .xls), or text (.txt, .csv, .dat)


    # module_input = f1

    crule = module_input["crule"]
    arule = module_input["arule"] #"/home/christos/Desktop/SCiO_Projects/REGALE/landmapr/LITAP/inst/extdata/arule.dbf"

    with open(module_input["input_json"], 'r') as file_:
          module_input = json.load(file_) 

    # folder = sys.argv[1] #"../../python_outputs/"
    use_spark = module_input["hyperparameters"]["use_spark"]
    end = module_input["hyperparameters"]["end"]
    n_remove = module_input["hyperparameters"]["n_remove"]
    procedure = module_input["hyperparameters"]["procedure"]
    zone = module_input["hyperparameters"]["zone"]
    out_format = module_input["hyperparameters"]["out_format"]
    clean = module_input["hyperparameters"]["clean"]
    verbose = module_input["hyperparameters"]["verbose"]
    quiet = module_input["hyperparameters"]["quiet"]
    log = module_input["hyperparameters"]["log"]
    nrow = module_input["hyperparameters"]["nrow"]
    ncol = module_input["hyperparameters"]["ncol"]


    resume = ""

    # TMP_DIR = "/tmp"
    TMP_DIR = "../data/"

    # output_backup_folder = TMP_DIR + "/python_outputs/backup/"
    # output_stats_folder = TMP_DIR + "/python_outputs/facet/"

    out_directory = TMP_DIR + "test_dem_outputs/"

    #%% 
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


    #%% Load previous Data from flow_mapper and form_mapper

    # get fill dem from flow_mappper
    db = get_previous(out_directory,step="fill",where="flow")
    db = dplyr.select(db, f.seqno, f.row, f.col, f.elev, f.drec, f.upslope, f.fill_shed, f.local_shed)
    db = add_buffer(db)

    # get form dem from flow_mappper
    weti = get_previous(out_directory,step="weti",where="form")
    weti = weti.drop(columns=dplyr.tidyselect.any_of(weti,["seqno_buffer", "drec_buffer"]))
    weti = weti.rename(columns={"qweti1" : "qweti", "qarea1" : "qarea", "lnqarea1" : "lnqarea"})
    weti = add_buffer(weti)

    relief = get_previous(out_directory, step = "relief", where = "form")
    relief = relief.drop(columns=dplyr.tidyselect.any_of(relief,["seqno_buffer", "drec_buffer"]))
    relief = add_buffer(relief)

    os.system("mkdir " + out_directory + "facet/")


    #%% Load/Derive the arule and crule information from input data or the execution data
    arule=None
    # Get Rules ---------------------------------------------------------------
    if arule is None:
        arule = arule_derive(weti, relief, n_remove = n_remove)
    else:
        # arule = "/home/christos/Desktop/SCiO_Projects/REGALE/regale/data/arule.csv"
        afile = arule
        # arule = load_extra(arule,type="arule") instead just load the csv
        arule = pd.read_csv(afile)
        arule.columns = [x.lower() for x in arule.columns]

        
    arule = format_rule(arule,"arule")

    # cfile = crule
    # folder = "~/Desktop/SCiO_Projects/REGALE/regale-ryax-modules/facet_mapper_version_2/"
    # cfile = folder + "data/crule.csv"
    crule = pd.read_csv(crule)
    crule.columns = [x.lower() for x in crule.columns]
    crule = format_rule(crule,"crule")


    check_rules(arule, crule)

    crule = prep_rule(crule,type_="crule")

    #%% Get Zones (if applicable)
    if procedure=="bc_pem":
        raise("TODO")
    else:
        weti["zone"] = 0
        arule["zone"] = 0
        crule["zone"] = 0

# #%%
# import findspark
# findspark.init()
# t = findspark.find()
# t
# #%%
# import os
# import sys
# from pyspark.sql import SparkSession
# import pyspark.sql.types
# from pyspark.sql.functions import udf, struct
# from pyspark.conf import SparkConf

# spark = SparkSession.builder.appName('pySparkSetup').getOrCreate()

# data_values = [('Apple',3),('Banana',6),('Orange', 9)]
# column_name = ['Name', 'Count']
# df = spark.createDataFrame(data_values).toDF(*column_name)
# df.show()

# os.environ['SPARK_HOME'] = r"C:\Users\Christos\anaconda3\envs\spark_env\Lib\site-packages\pyspark"
# os.environ['PYSPARK_PYTHON'] = r"C:\Users\Christos\anaconda3\envs\spark_env\python.exe"
# # os.environ['PYSPARK_DRIVER_PYTHON'] = r"C:\Users\Christos\anaconda3\envs\regale_spark\Lib\site-packages\pyspark"
# # os.environ['JAVA_HOME'] = r"C:\Program Files\Java\jdk-11"
# # os.environ['SPARK_LOCAL_IP'] = "localhost"

# # = C:\Program Files\Java\jdk1.8.0_201

# # Import PySpark


# #Create SparkSession
# spark = SparkSession.builder.appName('SparkByExamples.com').getOrCreate()

# # Data
# data = [("Java", "20000"), ("Python", "100000"), ("Scala", "3000")]

# # Columns
# columns = ["language","users_count"]

# # Create DataFrame
# df = spark.createDataFrame(data).toDF(*columns)

# # Print DataFrame
# df.show()


    #%% Facets

    if resume=="" or resume=="attributes":
        # Get attributes
        attr =  get_attr(weti,relief)
        
        # Create holder data
        fuzzattr = attr[['seqno', 'new_asp']]
        
        ##--------------CALL TO THE FUNCTION THAT HAS SPARK CODE--------------##
        ##Get fuzzy attributes
        ## Change use_spark argument in order to enable and disable spark and compare execution times
        fuzzattr = lsm_fuza(attr = attr, arule = arule, procedure = procedure, use_spark=use_spark)
        ##--------------------------------------------------------------------##
        
        save_output2(data=fuzzattr, name="fuza", locs=out_directory, out_format=out_format, where = "facet", add_db=db[["seqno", "buffer", "row", "col"]])


    #%%
    # Calculating Classes

    # get_previous but not sure if needed because will not start from this point on 
    # most probably but will be a continuous flow of the script
    first_element_not_buffer = (db["buffer"]==False).idxmax()
    fuzzattr_C = lsm_fuzc(fuzzattr, crule,first_element_not_buffer)

    save_output2(data=fuzzattr_C, name="fuzc", locs=out_directory, out_format=out_format, where = "facet", add_db=db[["seqno", "buffer", "row", "col"]])
    
    #%%
    
    classes = np.unique(crule["f_name"]).tolist()
    df_for_tif = fuzzattr_C.loc[(fuzzattr_C[classes].sum(axis=1) != 0)]
    
    # nrow = 327
    # ncol = 274
    # nrow = 339
    # ncol = 265
    
    first_option_array = np.reshape(df_for_tif["max_facet"].to_numpy(), (ncol,nrow))

    labels_transition_dictionary = {}
    for i,uniq in enumerate(np.unique(first_option_array)):
        labels_transition_dictionary[uniq] = i
        first_option_array = np.where(first_option_array==uniq,labels_transition_dictionary[uniq],first_option_array)

    first_option_array = first_option_array.astype(np.int16)    

    # %% CREATE A FIGURE AND SAVE A TIF FILE
    # plt.figure(figsize=(20,20))
    # plt.imshow(first_option_array)

    output_file = out_directory + "classification_DEM_map.tif"
    imwrite(output_file, first_option_array)

    return {'output_file' : output_file}


#%%

# ../../final/python_epirus_2_example_outputs/facet_epirus_2_input_json.json
f1 = {"input_json" : ".../data/facet_test_dem_input_json.json",
  "arule" :"../data/arule.csv" ,
  "crule" :"../data/crule.csv"
  }

start_time = time.time()
print(start_time)
t = handle(f1)  
end_time = time.time()
print(end_time)
print("Total time:", str(end_time-start_time))


