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
import json
import requests
import time

#import from other files
from LITAP_utils import *
from LITAP_load import *
from LITAP_read_write import *
from LITAP_functions import *

import boto3
import logging
from botocore.exceptions import ClientError

def calc_relief(relz):
    r = copy.deepcopy(relz)
    r["min_elev"] = np.nanmin(r["elev"])
    r["max_elev"] = np.nanmax(r["elev"])
    r["elev_range"] = r["max_elev"] - r["min_elev"]
    # r["max_elev_shed"] = r["elev"].groupby(r['fill_shed'],group_keys=False).transform(max_na)
    r["max_elev_shed"] = r["elev"].groupby(r['fill_shed'],group_keys=False).transform(np.nanmax)
    
    r["zpit2peak"] = r["z2pit"] + r["z2peak"]
    r["zcr2st"] = r["z2cr"] + r["z2st"]
    r["zcr2pit"] = r["z2cr"] + r["z2pit"]
    r["z2top"] = r["max_elev_shed"] - r["elev"]
    r["ztop2pit"] = r["max_elev_shed"] - r["pit_elev"]
    r["ncr2st"] = r["n2cr"] + r["n2st"]
    r["pctz2top"] = datar.all.trunc(100 - (r["z2top"] / r["ztop2pit"]) * 100)
    r["pctz2st"] = datar.all.trunc((r["z2st"] / r["zcr2st"]) * 100)
    r["pctz2pit"] = datar.all.trunc((r["z2pit"] / r["zpit2peak"]) * 100)
    r["pctn2st"] = datar.all.trunc((r["n2st"] / r["ncr2st"]) * 100)
    r["pmin2max"] = datar.all.trunc(((r["elev"] - r["min_elev"]) / r["elev_range"]) * 100)
    
    
    r.loc[r["ztop2pit"]==0,"pctz2top"] = 0
    r.loc[r["zcr2st"]==0,"pctz2st"] = 0
    r.loc[r["ncr2st"]==0,"pctn2st"] = 0
    r.loc[r["elev_range"]==0,"pmin2max"] = 0
    
    return r

def calc_relz(db, idb, str_val = 10000, ridge_val = 10000, pond = None, verbose = False, path_to_use=None):
    #-------------SET TMP FODLER PATH---------------#
    tmp_folder = path_to_use
    #-------------SET TMP FODLER PATH---------------#
    if verbose:
        print("  Calculating streams")

    streams = copy.deepcopy(db)
    streams["shedno"] = db["fill_shed"]

    json_file = {
          "db" : streams.replace(np.nan,"null").to_dict(orient="index"),
          "str_val" : str_val
          }
    
    with open(tmp_folder + "python.json", 'w') as f:
        json.dump(json_file, f)
        
    # AWS TEST
    path_to_file_for_upload = tmp_folder + "python.json"
    # target_bucket = "r-lambdas-dummy"
    target_bucket = "r-lambdas-dummy"
    
    string = path_to_file_for_upload.split("/")
    object_name = string[-1]
    
    # Upload the file
    s3_client = boto3.client('s3',
                            aws_access_key_id="VALID_AWS_KEY", 
                            aws_secret_access_key="VALID_AWS_SECRET_KEY", 
                            region_name="eu-central-1"
                            )
    try:
        response = s3_client.upload_file(path_to_file_for_upload, target_bucket, object_name)
    except ClientError as e:
        logging.error(e)
    
    
    json_file = {"body": "https:/r-lambdas-dummy.s3.eu-central-1.amazonaws.com/python.json"}
    #JSON FILE IS NOT USED AT ALL.. EVERYTHING IS HARDCODED
    response = requests.post("https://lambda.regale.form-03-calc-relz.calc-stream3.scio.services/", json = json.dumps(json_file))
    print(1, response)
    
    time.sleep(30)
    
    response = requests.get("https://r-lambdas-dummy.s3.eu-central-1.amazonaws.com/output.json")
    print(response)
    streams = pd.DataFrame(response.json())
    
    #---------------------------------------------------------------------------------------------------------------
    if verbose:
        print("  Calculating stream to pits")
        
    str2pits = copy.deepcopy(db)
    str2pits["shedno"] = db["local_shed"]
    
    json_file = {
         "db" : str2pits.replace(np.nan,"null").to_dict(orient="index"),
         "str_val" : ridge_val,
         "pond" : pond.to_dict(orient="index")
         }
    
    with open(tmp_folder + "python_A.json", 'w') as f:
        json.dump(json_file, f)
    
    # AWS TEST
    path_to_file_for_upload = tmp_folder + "python_A.json"
    # target_bucket = "r-lambdas-dummy"
    
    string = path_to_file_for_upload.split("/")
    object_name = string[-1]
    
    try:
        response = s3_client.upload_file(path_to_file_for_upload, target_bucket, object_name)
    except ClientError as e:
        logging.error(e)        
    
    json_file = {"body": "https://r-lambdas-dummy.s3.eu-central-1.amazonaws.com/python_A.json"}
    
    response = requests.post("https://lambda.regale.form-03-calc-relz.calc-pit3.scio.services/", json = json_file)
    print(2, response)
    time.sleep(30)
    
    response = requests.get("https://r-lambdas-dummy.s3.eu-central-1.amazonaws.com/output_A.json")
    
    str2pits = pd.DataFrame(response.json())   
    # return str2pits
    #---------------------------------------------------------------------------------------------------------------
    
    
    if verbose:
        print("  Calculating ridges") 
   
    json_file = {
         "db" : idb.replace(np.nan,"null").to_dict(orient="index"),
         "str_val" : ridge_val
         }
    
    with open(tmp_folder + "python.json", 'w') as f:
        json.dump(json_file, f)
        
        
    # AWS TEST
    path_to_file_for_upload = tmp_folder + "python.json"
    # target_bucket = "r-lambdas-dummy"
     
    string = path_to_file_for_upload.split("/")
    object_name = string[-1]
     
    try:
        response = s3_client.upload_file(path_to_file_for_upload, target_bucket, object_name)
    except ClientError as e:
        logging.error(e)        
     
    json_file = {"body": "https://r-lambdas-dummy.s3.eu-central-1.amazonaws.com/python.json"}
        
    response = requests.post("https://lambda.regale.form-03-calc-relz.calc-stream3.scio.services/", json = json_file)
    print(3, response)
    time.sleep(30)
    
    response = requests.get("https://r-lambdas-dummy.s3.eu-central-1.amazonaws.com/output.json")
    
    ridges = pd.DataFrame(response.json())
    
    rename_columns_dict = {
        "str_row" : "cr_row",
        "str_col" : "cr_col",
        "str_elev" : "cr_elev",
        "z2st" : "z2cr",
        "n2st" : "n2cr"
        }
    ridges = ridges.rename(columns=rename_columns_dict)
    ridges["cr_elev"] = np.nanmax(db["elev"]) - ridges["cr_elev"]
    
    
    #---------------------------------------------------------------------------------------------------------------
    if verbose:
        print("  Calculating ridges to pits") 
   
    ridge2pits = copy.deepcopy(idb)
    
    json_file = {
         "db" : ridge2pits.replace(np.nan,"null").to_dict(orient="index"),
         "str_val" : ridge_val,
         "pond" : pd.DataFrame().to_dict(orient="index")
         }
    
    with open(tmp_folder + "python_A.json", 'w') as f:
        json.dump(json_file, f)
        
        
    # AWS TEST
    path_to_file_for_upload = tmp_folder + "python_A.json"
    # target_bucket = "r-lambdas-dummy"
     
    string = path_to_file_for_upload.split("/")
    object_name = string[-1]
     
    try:
        response = s3_client.upload_file(path_to_file_for_upload, target_bucket, object_name)
    except ClientError as e:
        logging.error(e)        
     
    json_file = {"body": "https://r-lambdas-dummy.s3.eu-central-1.amazonaws.com/python_A.json"}    
    
    response = requests.post("https://lambda.regale.form-03-calc-relz.calc-pit3.scio.services/", json = json_file)
    print(4, response)
    time.sleep(30)
    
    response = requests.get("https://r-lambdas-dummy.s3.eu-central-1.amazonaws.com/output_A.json")
    ridge2pits = pd.DataFrame(response.json())   
   
    rename_columns_dict = {
        "pit_seqno" : "peak_seqno",
        "pit_row" : "peak_row",
        "pit_col" : "peak_col",
        "pit_elev" : "peak_elev",
        "z2pit" : "z2peak",
        "n2pit" : "n2peak",
        }
    ridge2pits = ridge2pits.rename(columns=rename_columns_dict)
    ridge2pits["peak_elev"] = np.nanmax(db["elev"]) - ridge2pits["peak_elev"]
    
    
    all_ = dplyr.left_join(streams, str2pits, by="seqno")
    all_ = dplyr.left_join(all_, ridges, by="seqno")
    all_ = dplyr.left_join(all_, ridge2pits, by="seqno")
    
    temp = db[["seqno", "elev", "row", "col", "buffer", "fill_shed"]]
    all_ = dplyr.left_join(all_, temp, by="seqno")
    
    all_.loc[np.isnan(all_["elev"]),["z2st","z2cr"]] = np.nan
    
    #---------------------------------------------------------------------------------------------------------------
    if verbose:
        print("  Calculating ridges to pits") 
        
    # return all_

    all_ = calc_relief(all_)        
    
    return all_
