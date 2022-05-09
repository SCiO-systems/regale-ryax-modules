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



def calc_relief(relz):
    r = copy.deepcopy(relz)
    r["min_elev"] = np.nanmin(r["elev"])
    r["max_elev"] = np.nanmax(r["elev"])
    r["elev_range"] = r["max_elev"] - r["min_elev"]
    r["max_elev_shed"] = r["elev"].groupby(r['fill_shed'],group_keys=False).transform(max_na)
    
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

def calc_relz(db, idb, str_val = 10000, ridge_val = 10000, pond = None, verbose = False):
    if verbose:
        print("  Calculating streams")

    streams = copy.deepcopy(db)
    streams["shedno"] = db["fill_shed"]

    json_file = {
          "db" : streams.replace(np.nan,"null").to_dict(orient="index"),
          "str_val" : str_val
          }
    
    with open("/home/christos/Desktop/SCiO_Projects/REGALE/regale/code/python.json", 'w') as f:
        json.dump(json_file, f)
        
    response = requests.post("https://lambda.regale.form-03-calc-relz.calc-stream3.scio.services/", json = json_file)
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
    
    with open("/home/christos/Desktop/SCiO_Projects/REGALE/regale/code/python_A.json", 'w') as f:
        json.dump(json_file, f)
        
    response = requests.post("https://lambda.regale.form-03-calc-relz.calc-pit3.scio.services/", json = json_file)
    
    str2pits = pd.DataFrame(response.json())   
    
    #---------------------------------------------------------------------------------------------------------------
    
    
    if verbose:
        print("  Calculating ridges") 
   
    json_file = {
         "db" : idb.replace(np.nan,"null").to_dict(orient="index"),
         "str_val" : ridge_val
         }
    
    with open("/home/christos/Desktop/SCiO_Projects/REGALE/regale/code/python.json", 'w') as f:
        json.dump(json_file, f)
        
        
    response = requests.post("https://lambda.regale.form-03-calc-relz.calc-stream3.scio.services/", json = json_file)
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
    
    with open("/home/christos/Desktop/SCiO_Projects/REGALE/regale/code/python_A.json", 'w') as f:
        json.dump(json_file, f)
        
    response = requests.post("https://lambda.regale.form-03-calc-relz.calc-pit3.scio.services/", json = json_file)
    
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
    
    all_ = calc_relief(all_)        
    
    return all_


# def calc_stream3(db, str_val = 10000, verbose = True):

#     l = [None]
#     str_row = pd.Series(len(db)*l)
#     str_col = pd.Series(len(db)*l)
#     str_elev = pd.Series(len(db)*l)
#     n2st = pd.Series(len(db)*l)
#     l = [0]
#     a_rep = pd.Series(len(db)*l)
#     b_rep = pd.Series(len(db)*l)
#     z2st = pd.Series(len(db)*l)
    
#     seqno_order = dplyr.select(db, f.seqno, f.row, f.col, f.elev, f. upslope, f.shedno)
#     seqno_order = dplyr.filter(seqno_order, ~pd.isna(seqno_order["elev"]))
#     seqno_order = dplyr.arrange(seqno_order,f.shedno , dplyr.desc(seqno_order["elev"]),f.upslope)
#     seqno_order = seqno_order["seqno"]
    
#     db_temp = dplyr.arrange(db, db["seqno"])
#     db_upslope = db_temp["upslope"]
#     db_row = db_temp["row"]
#     db_col = db_temp["col"]
#     db_elev = db_temp["elev"]
#     db_drec = db_temp["drec"]
    
#     n_cells = 250
#     seqno_track_groups = list(range(n_cells, len(seqno_order) + n_cells - 1, n_cells))
    
#     for idx,a in enumerate(seqno_track_groups):   
#         print(a)
#         # Remove extra cells on last batch
#         if a==seqno_track_groups[-1]:
#             seqno_sub = seqno_order[a-n_cells:seqno_track_groups[len(seqno_track_groups)-1]]
#         else:
#             seqno_sub = seqno_order[a-n_cells:a]
        
#         track_group = trace_flow_all(cells = seqno_sub, drec = db_drec)
#         # return track_group
            
#         for b,i in enumerate(seqno_sub):

                  

#             track = track_group[:,b].astype(int)
            
#             track,index = np.unique(track,return_index=True)
#             track = track[index].astype(int)
            
#             # Get first cell already visited     
#             visited = list(datar.base.which(z2st[track-1]>0))

#             if visited:
#                 visited = visited[0]
#             else:
#                 visited = 99999
            
#             # Get first channel cell
            
#             channel = list(datar.base.which(db_upslope.iloc[track-1] >= str_val))
#             if channel:
#                 channel = channel[0]
#             else:
#                 channel = 99999
      
#             # Get final cell          
#             end = track[min(visited, channel, len(track)-1)]
            
#             num_dn = datar.base.which(track==end)[0]

#             # Get new track
#             track = track[0:num_dn+1]
            
#             # Update channel and visited
#             visited = list(datar.base.which(z2st[track-1]>0))
#             channel = list(datar.base.which(db_upslope.iloc[track-1] >= str_val))
        
#             # If we go to the end (i.e. no visited, no channels)
#             # OR have a channel (and/or visited)
#             if not visited and not channel or channel:
#                 str_row[track-1] = db_row[end-1]
#                 str_col[track-1] = db_col[end-1]
#                 str_elev[track-1] = db_elev[end-1]
                
#                 a_rep[track-1] = a
#                 b_rep[track-1] = b + 1
                
#                 z2st[track-1] = db_elev.iloc[track-1] - db_elev.iloc[end-1]
#                 n2st[track-1] = list(range(num_dn,-1,-1))[:len(track)]
                
#                 # print("track: ", track)
#                 # print("end: ", end)
#                 # print("list: ", list(range(num_dn,-1,-1))[:len(track)])
#                 # print(not visited and not channel or channel)
#                 # print(z2st[track-1])
#                 # print("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")
                
#             else:
#                 # print("YEEEEEEEEEEEEEEEEEEEEEEES")
#                 str_row[track-1] = str_row[end-1]
#                 str_col[track-1] = str_col[end-1]
#                 str_elev[track-1] = str_elev[end-1]
                
#                 a_rep[track-1] = a
#                 b_rep[track-1] = b + 1
                    
#                 z2st[track-1] = db_elev.iloc[track-1] - str_elev.iloc[end-1]
#                 try:
#                     n2st[track-1] = list(range(num_dn + n2st.values[end-1],-1,-1))[:len(track)]
                    
                    
#                     # if idx==3 and len(track)>1:
#                     #     print("+++++++++++++++++++++++++++++++++")
#                     #     print(n2st)
#                     #     print(track)
#                     #     print(num_dn)
#                     #     print(n2st.values[end-1])
#                     #     print(list(range(num_dn + n2st.values[end-1] -1,-1,-1))[:len(track)])
                        
                        
#                 except Exception as e:
#                     print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
#                     print(visited)
#                     print(channel)
                    
#                     print("track: ", track)
#                     print("end: ", end)
#                     print("list: ", list(range(num_dn,-1,-1))[:len(track)])
#                     print(not visited and not channel or channel)
#                     print(z2st[track-1])

#                     raise(e)
#                     # print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
#                     # print(db_elev.iloc[track-1])
#                     # print(str_elev.iloc[end])
                
#             # if a==1000 and b==217:    
#             if a==1750 and b==10:
#                 print(visited)
#                 print(channel)
#                 # print()
#                 # print()
#                 # print()
                
#                 print("track: ", track)
#                 print("end: ", end)
#                 print("list: ", list(range(num_dn,-1,-1))[:len(track)])
#                 print(not visited and not channel or channel)
#                 print(z2st[track-1])
#                 return         
#             # if b==1:
#             #     break  
#         print("-------------------------------------------------")
#         # if a==250:
#         #     break  
        

#     z2st[z2st<0] = 0
    
#     df_dict = {
#         "seqno" : db_temp["seqno"],
#         "str_row" : str_row,
#         "str_col" : str_col,
#         "str_elev" : str_elev,
#         "z2st" : z2st,
#         "n2st" : n2st,
#         "a_rep" : a_rep,
#         "b_rep" : b_rep
#         }
    
#     return pd.DataFrame.from_dict(df_dict)



# def calc_pit3(db,pond=pd.DataFrame(),verbose=False):
#     if  len(pond.index)==0 or pond.empty:
#         # temp <- db %>%
#         # dplyr::filter(ddir == 5) %>%
#         # dplyr::select(shedno, pit_seqno = seqno, pit_row = row, pit_col = col, pit_elev = elev)
#         temp = dplyr.filter(db,f.ddir==5)
#         temp = dplyr.select(temp, f.shedno, pit_seqno = f.seqno, pit_row = f.row, pit_col = f.col, pit_elev = f.elev)
#     else:
#         temp = dplyr.select(pond, f.shedno, f.pit_seqno, f.pit_row, f.pit_col, f.pit_elev)
    
    
#     temp = dplyr.full_join(dplyr.select(db, f.seqno, f.shedno, f.elev), temp, by="shedno")
#     temp["z2pit"] = temp["elev"] - temp["pit_elev"]
#     temp["z2pit"].clip(lower=0)
#     temp["n2pit"] = pd.NA
#     temp.sort_values(by = 'seqno',inplace=True,ignore_index=True)
#     temp = temp.dropna(subset=["seqno"])
    
#     seqno_order = dplyr.arrange(db, db["upslope"], dplyr.desc(db["elev"]))
    
#     seqno_order = seqno_order.loc[~seqno_order["shedno"].isnull()]
#     seqno_order = seqno_order.loc[~seqno_order["elev"].isnull()]
    
#     seqno_order = seqno_order["seqno"]
    
    
#     n2pit = temp["n2pit"]
    
#     db_drec = db["drec"]
#     n_cells = 250
    
#     seqno_track_groups = list(range(n_cells, len(seqno_order) + n_cells, n_cells))
    
#     for a in seqno_track_groups:
            
#         # Remove extra cells on last batch
#         if a==seqno_track_groups[-1]:
#             seqno_sub = seqno_order.iloc[a-n_cells:]
#         else:
#             seqno_sub = seqno_order.iloc[a-n_cells:a]
        
#         track_group = trace_flow_all(cells=seqno_sub, drec=db_drec)
        
#         for b,i in enumerate(seqno_sub):
    
#             if pd.isna(n2pit.iloc[i]):
#                 track = track_group[:,b]
#                 track,index = np.unique(track,return_index=True)
#                 track = track[index].astype(int)
#                 n2pit.loc[track] = list(range(len(track)-1,-1,-1))
    
    
#     n2pit = n2pit.drop(index=0)
#     n2pit = pd.concat([n2pit,pd.DataFrame([pd.NA])],ignore_index=True)
#     temp["n2pit"] = n2pit.values
#     temp = dplyr.select(temp,~f.shedno,~f.elev)

#     return temp