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


def trace_wetness(db_n_sub, db_c):
    
    # Get cells with only lower cells
    temp = db_c[np.logical_and(db_c["in_t"] == 0,db_c["out_t"]>0)]
    
    
    # Get the lower cells that these upper cells flow into, add effects of this link,
    # plus all previous linkages from this upper cell
    temp_n = db_n_sub.loc[db_n_sub["seqno_n"].isin(temp["seqno"])]
    
    #next 3 command do the left join
    temp_n = temp_n.merge(temp,how='left',left_on="seqno_n",right_on="seqno")
    temp_n = temp_n.drop(columns=['seqno_y'])
    temp_n = datar.all.rename(temp_n, seqno = "seqno_x")
    
    temp_n["qc1"] = dplyr.if_else(temp_n["status"]=="higher", temp_n["qarea1"]/temp_n["sumtanbl"], temp_n["qarea1"]/0.0001)
    temp_n["qc2"] = dplyr.if_else(temp_n["status"]=="higher", temp_n["qarea2"]/temp_n["sumtanbl"], temp_n["qarea2"]/0.0001)
    
    # -1 because calc from lower cell's perspective
    temp_n["new_qa1"] = dplyr.if_else(temp_n["status"]=="higher", temp_n["qc1"]*(-1)*temp_n["tan2_f"], temp_n["qarea1"])
    temp_n["new_qa2"] = dplyr.if_else(temp_n["status"]=="higher", temp_n["qc2"]*(-1)*temp_n["tan2_f"], temp_n["qarea2"])
    
    temp_n.loc[temp_n["status"]!="higher","elev_diff"] = 0
    temp_n.loc[temp_n["status"]!="higher","cell_t"] = 0
    
    
    
    temp_qc = dplyr.group_by(temp_n, f.seqno_n)
    temp_qc = dplyr.arrange(temp_qc, f.order_n)
    # get last qc for each higher cell
    temp_qc = dplyr.summarize(temp_qc, 
                           qc1 = dplyr.last(f.qc1),
                           qc2 = dplyr.last(f.qc2))
    temp_qc = dplyr.arrange(temp_qc, f.seqno_n)
    
    
    temp_n = dplyr.group_by(temp_n, f.seqno)
    temp_n = dplyr.summarize(temp_n,
                           qarea1 = datar.all.sum_(f.new_qa1),
                           qarea2 = datar.all.sum_(f.new_qa2),
                           this_round = datar.all.sum_(f.status.isin(["higher", "higher_drec"])),
                           d = datar.all.sum_(f.diag) + datar.all.sum_(f.count_d),
                           o = datar.all.sum_(~f.diag) + datar.all.sum_(f.count_o),
                           t = f.d + f.o,
                           elev_sum = datar.all.sum_(f.elev_diff) + datar.all.sum_(f.elev_sum))

    
    temp_seqno = temp["seqno"]
    temp_n_seqno = temp_n["seqno"]
    temp_n_t = temp_n["t"]
    temp_n_d = temp_n["d"]
    temp_n_o = temp_n["o"]
    temp_n_elev_sum = temp_n["elev_sum"]
    temp_n_qarea1 = temp_n["qarea1"]
    temp_n_qarea2 = temp_n["qarea2"]
    temp_n_this_round = temp_n["this_round"]
    
    temp_qc_seqno1 = temp_qc["seqno_n"]
    temp_qc_seqno2 = temp_qc["seqno_n"]
    temp_qc_qc1 = temp_qc["qc1"]
    temp_qc_qc2 = temp_qc["qc2"]
    
    
    # Add these cell totals to all the previous totals already calculated for this lower cell
    seqno = db_c["seqno"]
    cell_t = db_c["cell_t"]
    count_d = db_c["count_d"]
    count_o = db_c["count_o"]
    elev_sum = db_c["elev_sum"]
    qarea1 = db_c["qarea1"]
    qarea2 = db_c["qarea2"]
    qc1 = db_c["qc1"]
    qc2 = db_c["qc2"]
    
    in_t = db_c['in_t']
    
    i = seqno.isin(temp_n_seqno)
    
    
    cell_t.loc[i] += temp_n_t.values
    count_d.loc[i] += temp_n_d.values
    count_o.loc[i] += temp_n_o.values
    elev_sum.loc[i] += temp_n_elev_sum.values
    qarea1.loc[i] += temp_n_qarea1.values
    qarea2.loc[i] += temp_n_qarea2.values
    
    qc1.loc[seqno.isin(temp_qc_seqno1)] = temp_qc_qc1.values
    qc2.loc[seqno.isin(temp_qc_seqno2)] = temp_qc_qc2.values
    
    # Resolve the linkages (i.e. one less upper cell flowing into these lower cells)
    in_t.loc[seqno.isin(temp_seqno)] -= 1
    in_t.loc[i] -= temp_n_this_round.values
    
    dataframe_dict = {
        "seqno" : seqno,
        "drec" : db_c["drec"],
        "cell_t" : cell_t,
        "count_d" : count_d,
        "count_o" : count_o,
        "elev_sum" : elev_sum,
        "qarea1" : qarea1,
        "qarea2" : qarea2,
        "qc1" : qc1,
        "qc2" : qc2,
        "in_t" : in_t,
        "out_t" : db_c["out_t"],
        "sumtanbl" : db_c["sumtanbl"]
        }
    
    return pd.DataFrame(data=dataframe_dict)


def calc_weti(db, grid = 5, verbose = False):
    l1 = grid * 0.5    # orthogonal
    l2 = grid * 0.354  # diagonal  (hypothenus = sqrt(0.5^2 + 0.5^2) / 2)
    l_sqr = grid * grid  # let's use cell area
    orthogonal = grid
    diagonal = grid * np.sqrt(2)
    qarea1 = 1
    qarea2 = l_sqr
    
    max_s = db.shape[0]
    
    if(verbose):
        print("  Setting up cell flow")
    
    
    
    db_n = dplyr.arrange(db, dplyr.desc(f.elev), f.upslope, f.seqno)
    db_n["order"] = list(range(1,db_n.shape[0]+1))
    
    db_n = dplyr.select(db_n, f.seqno, f.elev, f.drec, f.order)
    
    db_n = dplyr.arrange(db_n, f.seqno)
    
    db_n = nb_values(db=db_n, max_cols=max(db["col"]), col=["elev", "seqno", "order"], format_="wide")

    n1=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n1")
    n2=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n2")
    n3=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n3")
    n4=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n4")
    n5=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n5")
    n6=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n6")
    n7=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n7")
    n8=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n8")
    n9=["seqno","elev","drec","order"] + dplyr.tidyselect.contains(db_n, "n9")
    
    # MAYBE THIS CAN BE PARALLELIZED
    #---------------------------------
    db_n = tidyr.nest(db_n, n1=n1,n2=n2,n3=n3,n4=n4,n5=n5,n6=n6,n7=n7,n8=n8,n9=n9)
                            
    db_n = tidyr.pivot_longer(db_n, cols = dplyr.tidyselect.everything())
    
    #this is the big and complex dplyr.mutate
    
    # add column n in the first level of nested data
    db_n["n"] = list(range(1,10))
    
    for i in range(9):
        # remove the number from the columns name inside the second level of nested data
        name_dict = {}
        for col in db_n["value"].iloc[i].columns:
            if "_n" in col:
                name_dict[col] = col[:-1]
        db_n["value"].iloc[i].rename(columns=name_dict,inplace=True)
        
        # add the of the first level as a columns in the second level
        db_n["value"].iloc[i]["n"] = i+1
        
        # deternine if the data of the second level belong to diagonal elements or not
        db_n["value"].iloc[i]["diag"] = i+1 in [1,3,7,9]
        
        # create deltax value based on the True/False diag
        db_n["value"].iloc[i]["deltax"] = dplyr.if_else(db_n["value"].iloc[i]["diag"], diagonal, orthogonal)
        # create ql value based on the True/False diag
        db_n["value"].iloc[i]["ql"] = dplyr.if_else(db_n["value"].iloc[i]["diag"], l2, l1)
        
        db_n["value"].iloc[i]["status"] = dplyr.case_when(db_n["value"].iloc[i]["elev_n"]>db_n["value"].iloc[i]["elev"],"higher",
                                                          db_n["value"].iloc[i]["elev_n"]<db_n["value"].iloc[i]["elev"],"lower",
                                                          True,"no_flow")
        
        db_n["value"].iloc[i]["status_drec"] = np.logical_and(db_n["value"].iloc[i]["drec"]==db_n["value"].iloc[i]["seqno_n"],
                                                              db_n["value"].iloc[i]["drec"] !=db_n["value"].iloc[i]["seqno"])
        
        db_n["value"].iloc[i]["elev_diff"] = db_n["value"].iloc[i]["elev_n"] - db_n["value"].iloc[i]["elev"] 
        
        db_n["value"].iloc[i]["tan_f"] = -db_n["value"].iloc[i]["elev_diff"] / db_n["value"].iloc[i]["deltax"]
        
        db_n["value"].iloc[i]["tan2_f"] = db_n["value"].iloc[i]["tan_f"] * db_n["value"].iloc[i]["ql"]
        
    #---------------------------------
    
    db_n = dplyr.select(db_n, f.value)
    db_n = tidyr.unnest(db_n, f.value)
    
    
    
    if(verbose):
        print("  Addressing special flow cases")
    
    # Which cells flow to drec (not by elev)
    db_drec = dplyr.group_by(db_n, f.seqno)
    db_drec = db_drec.groupby("seqno").filter(lambda x: all(x["status"]!="lower"))
    db_drec = db_drec.groupby("seqno").apply(lambda x: x[x["status_drec"]])
    
    
  
    # Which should be evaluated BEFORE their drain points?
    db_first = dplyr.filter(db_drec, f.order<f.order_n)

    db_n.loc[np.logical_and(db_n["seqno"].isin(db_first["seqno"]),db_n["status_drec"]),"status"] = "lower_drec"

    # print(len(db_n.loc[db_n["status"]=="higher"]))

    db_n["match"] = db_n["seqno"].apply(str) +  "_" + db_n["seqno_n"].apply(str)
    db_n["match"] = db_n["match"].apply(lambda x: x.split(".")[0])

    db_first["match"] = db_first["drec"].apply(str) +  "_" + db_first["seqno"].apply(str)
    db_first["match"] = db_first["match"].apply(lambda x: x.replace(".0",""))

    db_n.loc[db_n["match"].isin(db_first["match"]),"status"] = "higher_drec"

    db_after = db_drec.loc[~db_drec["seqno"].isin(db_first["seqno"]),"seqno"]
    
    # Only cells which have higher neighbours
    db_n_sub = dplyr.select(db_n, f.seqno, f.seqno_n, f.order, f.order_n, f.diag, f.elev_diff, f.status, f.tan2_f)
    db_n_sub = db_n_sub.replace(np.nan,pd.NA)
    db_n_sub = dplyr.filter(db_n_sub, ~pd.isna(db_n_sub["elev_diff"]),db_n_sub["status"].isin(["higher", "higher_drec"]))

    
    db_n = db_n.replace(np.nan,pd.NA)
    db_c = dplyr.filter(db_n, ~pd.isna(db_n["elev_diff"]))
    db_c = dplyr.group_by(db_c, f.seqno)
    db_c = dplyr.summarize(db_c, 
                           drec = f.drec[0],
                           in_t = datar.all.sum_(f.status.isin(["higher", "higher_drec"])),
                           out_t = datar.all.sum_(f.status.isin(["lower", "lower_drec"])),                       
                           qarea1 = qarea1,
                           qarea2 = qarea2,
                           cell_t = 0,
                           count_d = 0,
                           count_o = 0,
                           elev_sum = 0,
                           sumtanbl = datar.all.sum_(f.tan2_f[f.status == "lower"]),
                           qc1 = pd.NA,
                           qc2 = pd.NA)
    
    
    if(verbose):
        print("  Trace wetness")
        
    #THIS WHILE DOES 1 LOOP LESS THAN R CODE, IF PROBLEM HAPPENS CHECK IT OUT
    i = 1
    while any(np.logical_and(db_c["in_t"]>=0,db_c["out_t"]!=0)):
        db_c = trace_wetness(db_n_sub, db_c)
        # print(i)
        i += 1

    if(verbose):
        print("  Fix special cases")        

    db_flat = dplyr.filter(db_c, db_c["seqno"].isin(db_after))
    db_flat = dplyr.select(db_flat, f.drec, qarea_flat1 = f.qarea1, qarea_flat2 = f.qarea2)

    db_c_temp = dplyr.inner_join(dplyr.select(db_c, f.seqno, f.qarea1, f.qarea2), db_flat, by = {'seqno': 'drec'})
    db_c_temp = dplyr.group_by(db_c_temp, f.seqno, f.qarea1, f.qarea2)
    db_c_temp = dplyr.summarize(db_c_temp, 
                           qarea_flat1 = datar.all.sum_(f.qarea_flat1),
                           qarea_flat2 = datar.all.sum_(f.qarea_flat2))

    db_c_temp = dplyr.mutate(db_c_temp, qarea1 = f.qarea1 + f.qarea_flat1, qarea2 = f.qarea2 + f.qarea_flat2)

    db_c_temp = db_c_temp.drop(columns=["qarea_flat1","qarea_flat2"])


    db_c.loc[db_c["seqno"].isin(db_c_temp["seqno"]),"qarea1"] = db_c_temp["qarea1"].values
    db_c.loc[db_c["seqno"].isin(db_c_temp["seqno"]),"qarea2"] = db_c_temp["qarea2"].values

    db_c = db_c.replace(np.nan, pd.NA)
    db_c.loc[pd.isna(db_c["qc1"]),"qc1"] = db_c.loc[pd.isna(db_c["qc1"]), "qarea1"]/0.0001
    db_c.loc[pd.isna(db_c["qc2"]),"qc2"] = db_c.loc[pd.isna(db_c["qc2"]), "qarea2"]/0.0001
    #HERE THE BIG NUMBERS(GREATER THAN BILLIONS) ARE NOT THE SAME WITH R, MAYBE DUE TO OVERFLOW


    #R CODE COMMENT:
    ## Shouldn't this be all +1?    
    db_weti = copy.deepcopy(db_c)
    db_weti["qweti1"] = np.where(db_c["qc1"] > 1, np.log(db_c["qc1"].astype(float)), np.log(1 + db_c["qc1"].astype(float)))
    db_weti["qweti2"] = np.where(db_c["qc2"] > 1, np.log(db_c["qc2"].astype(float)), np.log(1 + db_c["qc2"].astype(float)))
    db_weti["qarea1"] = round(db_weti["qarea1"], 2)
    db_weti["qarea2"] = round(db_weti["qarea2"], 2)
    db_weti["qweti1"] = round(db_weti["qweti1"], 2)
    db_weti["qweti2"] = round(db_weti["qweti2"], 2)


    db = dplyr.full_join(db_weti, db, by = ['seqno', 'drec'])
    db = datar.all.arrange(db,f.seqno)
    
    
    return db