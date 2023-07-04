# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:54:22 2022

@author: yyzou
"""

import pandas as pd
import numpy as np

def readchr(chr_len_info,res,drop_chrom=None):
    chr_len = pd.read_table(chr_len_info,names=["chr","length"])
    if drop_chrom:
        chr_len = chr_len[~chr_len["chr"].astype('str').isin(drop_chrom)]
    
    chr_len["bins"] = round(chr_len["length"].astype("int")/int(res))
    chr_len["bin_num"] = [0] + [int(chr_len["bins"].loc[:i-1].sum()) for i in  chr_len.index[1:]]
    chr_len["chr"] = chr_len["chr"].astype("str")
    return chr_len



def fitbin(intera_mat,organ,cnv_info, chr_len, res):
    
    intera_mat["organ"] = organ

    intera_mat["count"]=intera_mat["count"].astype("int")
    intera_mat["start1"] =intera_mat["start1"].astype("int")
    intera_mat["start2"] =intera_mat["start2"].astype("int")
    chr_len["chr"] = chr_len['chr'].apply(lambda x: str(x))
    print(intera_mat["chrom2"].unique())
    intera_mat["chrom1_bin"] = intera_mat.apply(lambda x: int(x["start1"]/int(res))+int(chr_len["bin_num"].loc[chr_len["chr"]==str(x["chrom1"])])+1,axis=1)
    intera_mat["chrom2_bin"] = intera_mat.apply(lambda x: int(x["start2"]/int(res))+int(chr_len["bin_num"].loc[chr_len["chr"]==str(x["chrom2"])])+1,axis=1)
    
    ## region info
    cnv_region = pd.read_table(cnv_info,header=None,sep="\t")
    cnv_region["start_bin"] = cnv_region.apply(lambda x: int(int(x[1])/int(res))+int(chr_len["bin_num"].loc[chr_len["chr"]==str(x[0])])+1,axis=1)
    cnv_region["end_bin"] = cnv_region.apply(lambda x: int(int(x[2])/int(res))+int(chr_len["bin_num"].loc[chr_len["chr"]==str(x[0])])+1,axis=1)

    ### the interaction may presented as cnv_VS_costomer-region or costomer-region_VS_cnv or cnv_cnv 
    intera_mat["customer"] = intera_mat.apply(lambda x:x["chrom2_bin"] if x["chrom1_bin"] in range(cnv_region["start_bin"].loc[0],cnv_region["end_bin"].loc[0]+1) else x["chrom1_bin"],axis=1) 


    ### cumulativate count bin by bin
    
    cumulative_interaction_intensity = pd.DataFrame({"bin":range(1,int(chr_len["bins"].sum())),"count":0})
    cumulative_interaction_intensity["count"] =  cumulative_interaction_intensity.apply(lambda x: int(intera_mat["count"].loc[intera_mat["customer"]==x["bin"]].sum()) if x["bin"] in intera_mat["customer"].to_list() else 0, axis=1)
    cumulative_interaction_intensity["cumsum"] = np.cumsum(cumulative_interaction_intensity["count"])
    
    
    return cumulative_interaction_intensity

    
