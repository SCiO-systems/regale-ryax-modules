import sys
import os


os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(1, 'functions/')

#imports
# from dbfread import DBF
import pandas as pd
import numpy as np
# import copy
# import string
from datar import dplyr,tidyr
from datar.all import f
# from itertools import compress
# import math
import datar
import time

import pyspark.pandas as ps
import pyspark.sql.types
from pyspark.sql import SparkSession
from pyspark.sql.functions import udf, struct
from pyspark.conf import SparkConf



def prep_rule(rule, type_):
    if type_=='crule':
        r = rule.groupby(['f_name','zone']).apply(lambda x: x['attrwt']/np.sum(x['attrwt'])).reset_index()
        r = r.sort_values('level_2').reset_index(drop=True)
        rule['relwt'] = r['attrwt']
        rule = rule[['zone', 'f_name', 'fuzattr', 'relwt']]
    return rule


def get_attr(weti, relief):    
    arule_weti = ['prof', 'plan', 'slope_pct', 'aspect', 'qweti', 'qarea', 'lnqarea', 'new_asp']
    arule_relief = ['pctz2st', 'pctz2pit', 'z2pit']
    attr = weti[arule_weti + ['seqno','zone']]
    
    attr = dplyr.left_join(attr, relief[arule_relief+['seqno']],by='seqno')
    return attr

##----------------------------SPARK RELATED CODE START----------------------------
##--------------------------------------------------------------------------------

# NON-spark function for comparing execution times
def arule_models(model, x, b, b_low, b_hi, b1, b2, d):

    # Calculate fuzzy attributes for each cell
    def basic_model(x,b,d):
        return 1/(1+((x-b)/d)**2)

    if model==1:
        fuzz = basic_model(x, b, d)
    elif model==2:
        fuzz = dplyr.case_when(np.logical_and(x>b_low,x<b_hi),1,
                               x<=b_low,basic_model(x, b_low, d),
                               x>=b_hi,basic_model(x, b_hi, d))
        
    elif model==3:
        fuzz = dplyr.case_when(np.logical_and(x>b1,x<b2),1,
                               x<=b1,basic_model(x, b1, d),
                               x>=b2,basic_model(x, b2, d))
    elif model==4:
        fuzz = dplyr.if_else(x>b, 1, basic_model(x, b, d))
    elif model==5:
        fuzz = dplyr.if_else(x<b, 1, basic_model(x, b, d))
        
    return 100*fuzz

# SPARK function used with sparks' udf 
def arule_models_spark(row):
    
    model = row["model"]
    x = row["x"]
    b = row["b"]
    b_low = row["b_low"]
    b_hi = row["b_hi"]
    b1 = row["b1"]
    b2 = row["b2"]
    d = row["d"]

    # Calculate fuzzy attributes for each cell
    def basic_model(x,b,d):
        if d==0:
            d=1
        return 1/(1+((x-b)/d)**2)

    if model==1:
        fuzz = basic_model(x, b, d)
    elif model==2:
        fuzz = dplyr.case_when(np.logical_and(x>b_low,x<b_hi),1,
                               x<=b_low,basic_model(x, b_low, d),
                               x>=b_hi,basic_model(x, b_hi, d))
        
    elif model==3:
        fuzz = dplyr.case_when(np.logical_and(x>b1,x<b2),1,
                               x<=b1,basic_model(x, b1, d),
                               x>=b2,basic_model(x, b2, d))
    elif model==4:
        fuzz = dplyr.if_else(x>b, 1, basic_model(x, b, d))
    elif model==5:
        fuzz = dplyr.if_else(x<b, 1, basic_model(x, b, d))
        
    return float(100*fuzz)


# Function that calls the spark function to be applied on data
def lsm_fuza(attr, arule, procedure,use_spark=True):
    # Create holder data
    fuzzattr = attr[['seqno', 'new_asp']]
    
    ## condition for using spark in order to compare execution times
    if use_spark==True:

        #----------------------SPARK----------------------#
        # start a spark session
        spark = SparkSession.builder.appName("Spark_TEST").getOrCreate()
        #define as udf the function that must be applied to the data
        #the function is to be applied row-wise and return a float number as result
        arule_models_udf = udf(lambda row: arule_models_spark(row),returnType=pyspark.sql.types.FloatType())
        #start timing the execution
        start_time = time.time()
        ## PARALLEL
        ## this for-loop is necessary for the function must be applied 17 times on the data with different inputs each time
        ## at a later stage this for-loop can be absorbed in a spark workflow but for now the approach is this
        for a in range(len(arule)):
            print(a, arule['class_out'].iloc[a])
            #prepare the dataframe that are going to be fed to spark
            attr_zoned = dplyr.filter(attr, attr['zone']==arule['zone'].iloc[a])
            
            attr_zoned_ps = attr_zoned[arule['attr_in'].iloc[a]].to_frame("x")
            attr_zoned_ps["model"] = arule['model'].iloc[a]
            attr_zoned_ps["b"] = arule['b'].iloc[a]
            attr_zoned_ps["b_low"] = arule['b_low'].iloc[a]
            attr_zoned_ps["b_hi"] = arule['b_hi'].iloc[a]
            attr_zoned_ps["b1"] = arule['b1'].iloc[a]
            attr_zoned_ps["b2"] = arule['b2'].iloc[a]
            attr_zoned_ps["d"] = arule['d'].iloc[a]
             
            # create the spark dataframe 
            attr_zoned_ps = spark.createDataFrame(attr_zoned_ps)
            # call spark execution on the spark dataset created earlier
            attr_zoned_ps = attr_zoned_ps.withColumn(arule['class_out'].iloc[a], arule_models_udf(struct([attr_zoned_ps[x] for x in attr_zoned_ps.columns]))).collect()
            # transform the results to be compatible with local pandas dataframe manipulation
            attr_zoned_ps = spark.createDataFrame(attr_zoned_ps).toPandas()
            # extract the result column produced from the execution
            attr_zoned[arule['class_out'].iloc[a]] = attr_zoned_ps[arule['class_out'].iloc[a]]
            # select specific columns from the produced dataframe
            attr_zoned = attr_zoned[['seqno','zone',arule['class_out'].iloc[a]]]
            # save the selected columns to the final dataframe to be used later by the script
            fuzzattr.loc[attr_zoned['seqno']-1, attr_zoned.columns] = attr_zoned
        
        #stop the time            
        end_time = time.time()
        print("Time (seconds) to execute WITH SPARK: ", end_time - start_time) 
        #close the spark session to free resources
        spark.stop()
    else:
        start_time = time.time()
        ## NON-PARALLEL
        for a in range(len(arule)):
            print(a, arule['class_out'].iloc[a])
            attr_zoned = dplyr.filter(attr, attr['zone']==arule['zone'].iloc[a])
            attr_zoned[arule['class_out'].iloc[a]] = attr_zoned[arule['attr_in'].iloc[a]].apply(lambda row: arule_models(arule['model'].iloc[a], row, arule['b'].iloc[a], arule['b_low'].iloc[a], arule['b_hi'].iloc[a], arule['b1'].iloc[a], arule['b2'].iloc[a], arule['d'].iloc[a]).astype('float32'))
            attr_zoned[arule['class_out'].iloc[a]] = attr_zoned[arule['class_out'].iloc[a]].astype('float32')
            attr_zoned = attr_zoned[['seqno','zone',arule['class_out'].iloc[a]]]
            
            # fuzzattr[f$seqno, names(f)] <- f
            fuzzattr.loc[attr_zoned['seqno']-1, attr_zoned.columns] = attr_zoned
        
        end_time = time.time()
        print("Time (seconds) to execute WITHOUT SPARK: ", end_time - start_time) 
    
    if len([i for i in ['planar_d', 'planar_a'] if i in fuzzattr.columns])==2:
        fuzzattr['planar_2x'] = (fuzzattr['planar_d'] + fuzzattr['planar_a'])/2
        
    # For Second option
    if procedure=='bc_pem':
        raise('TODO in facet_01_lsm.R')
        
        
    return fuzzattr


##------------------------------------------------------------------------------
##----------------------------SPARK RELATED CODE END----------------------------


def fuzc_sum(fuzzattr, crule, first_element_not_buffer,use_spark=True):
    
    t = list(set(crule["fuzattr"]))

    #using d instead of f variable due to imports already using f value
    d = fuzzattr[["seqno", "zone"] + t]
    d = dplyr.arrange(d, f.zone, f.seqno)

    print("!checking zone for None values, example doesn't have any so in a real case here MAY be a BUG!")
    d = dplyr.filter(d, ~pd.isna(d["zone"]))

    seqnos = d["seqno"]

    # db_cleared = db.replace(pd.NA,-9999)

    d = tidyr.nest(d, data=~f.zone)
    d = dplyr.left_join(d, crule, by="zone")

    # attr_zoned[arule['class_out'].iloc[a]] = attr_zoned[arule['attr_in'].iloc[a]].apply(lambda row: arule_models(arule['model'].iloc[a], row, arule['b'].iloc[a], arule['b_low'].iloc[a], arule['b_hi'].iloc[a], arule['b1'].iloc[a], arule['b2'].iloc[a], arule['d'].iloc[a]).astype('float32'))
    
    # !!!!!!!TO BE PARALLELIZED!!!!!!!!!
    start_time = time.time()
    d["data"] = d.apply(lambda row: np.asarray(row["data"][row["fuzattr"]] * row["relwt"]),axis=1)#.astype('float32')
    end_time = time.time()
    print(end_time - start_time)

    r = d.groupby(["zone","f_name"]).apply(lambda x: np.sum(np.array(x["data"].tolist()),axis=0)).reset_index()
    r = r.rename(columns={0:"data"})
    r = pd.pivot_table(r,values="data",index="zone",columns="f_name").reset_index()

    k = tidyr.unnest(r, [x for x in r.columns if x!="zone"])
    
    k.index = np.arange(first_element_not_buffer, first_element_not_buffer + len(k))    
    d = pd.DataFrame(0, index=np.arange(len(fuzzattr)), columns=k.columns)

    d.iloc[k.index] = k
    d["seqno"] = seqnos

    return d


def fuzc_max(fuzc):
    #using D instead of f beracuse f is taken during imports
    d = fuzc[[x for x in fuzc.columns if x not in ["zone", "seqno"]]]
    n = d.columns.tolist()
    
    ## Slower approach
    # start_time = time.time()
    # #assign NA to all values per row except the 2 largest
    # d = d.apply(lambda row: row.nlargest(2),axis=1)
    # #get only the 2 max values and add as column the class from the index (next command will
    # # reverse the index but this is necessery for the command to be executed correctly)
    # d = d.apply(lambda row: row.loc[~np.isnan(row)].reset_index(), axis=1)
    # # assign the column of the classes as index again and transpose the dataframe
    # d = d.apply(lambda row: (row.set_index("index")).transpose())
    # #order the columns' position based on the value of the two max in order to formalize their position ("max_2nd_value", "max_value")
    # d = d.apply(lambda row: row[row.columns[row.loc[row.last_valid_index()].argsort()]])
    # # create the list combining the ordered max values and the respective columns
    # d = d.apply(lambda row: list(row.values[0]) + list(row.columns))
    # # create the final dataframe 
    # d = pd.DataFrame.from_dict(dict(zip(d.index, d.values)),orient="index",columns=["max_2nd_value", "max_value", "max_2nd_facet", "max_facet" ])
    # # end_time = time.time()
    # # print("1st approach total time: ", end_time-start_time)

    def max_2_classes_extractor(row):
        
        names = list(row.index[row.notnull()])
        values = list(row.dropna(axis=0, how='all'))
        
        zipped_lists = zip(values, names)
        sorted_pairs = sorted(zipped_lists)
    
        tuples = zip(*sorted_pairs)
        list1, list2 = [ list(tuple) for tuple in  tuples]
        
        
        return list1 + list2
    
    
    # start_time = time.time()
    d = d.apply(lambda row: max_2_classes_extractor(row.nlargest(2)),axis=1)
    d = pd.DataFrame.from_dict(dict(zip(d.index, d.values)),orient="index",columns=["max_2nd_value", "max_value", "max_2nd_facet", "max_facet" ])
    # end_time = time.time()

    # print("2nd approach total time: ", end_time-start_time)

    return dplyr.bind_cols(fuzc, d)

def lsm_fuzc(fuzzattr, crule,first_element_not_buffer,use_spark=True):
    print("YES")
    k = fuzc_sum(fuzzattr, crule,first_element_not_buffer,use_spark=True)
    # return k
    return fuzc_max(k)