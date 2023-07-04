# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:54:03 2022

@author: yyzou

Ussage: python zoom.py cnv.bed whole_genome_interaction.txt outname
"""



import pandas as pd
from sys import argv


def pre_cnv(cnv_file):
    cnv_region = pd.read_table(cnv_file,header=None,sep="\t")  ## bed format, [chrom, start, end, others]
    return cnv_region

def pre_mat(interaction_file,drop_chrom):
    whole_hic = pd.read_table(interaction_file)   ### "whole_genome_interaction.txt"
    if drop_chrom:     ### if not None
        whole_hic = whole_hic[(~whole_hic["chrom1"].astype('str').isin(drop_chrom)) & (~whole_hic["chrom2"].astype('str').isin(drop_chrom))]
    return whole_hic


def overlap(series,cnv):
    if len(cnv.shape)>1:
        for i in cnv.index:
            if str(cnv[0].loc[i])==str(series["chrom1"]):
                if (int(cnv[1].loc[i])-int(series["start1"]))*(int(cnv[2].loc[i])-int(series["end1"]))<=0:
                    return 1
            elif str(cnv[0].loc[i])==str(series["chrom2"]):
                if (int(cnv[1].loc[i])-int(series["start2"]))*(int(cnv[2].loc[i])-int(series["end2"]))<=0:
                    return 1
    else:
        if str(cnv[0])==str(series["chrom1"]):
            if (int(cnv[1])-int(series["start1"]))*(int(cnv[2])-int(series["end1"]))<=0:
                return 1
        elif str(cnv[0])==str(series["chrom2"]):
            if (int(cnv[1])-int(series["start2"]))*(int(cnv[2])-int(series["end2"]))<=0:
                return 1

def zoom(cnv_file, interaction_file,outname,drop_chrom=None):
    cnv_region = pre_cnv(cnv_file)
    whole_hic = pre_mat(interaction_file, drop_chrom)  ### data preparation

    whole_hic["overlap"] = whole_hic.apply(lambda x: overlap(x,cnv_region),axis=1)
    #### overlaped interactions
    target_interaction = whole_hic.loc[whole_hic["overlap"]==1]
    
    target_interaction = target_interaction.reset_index()
    target_interaction.to_csv(outname+"_target_region.xls",sep="\t",index=False)
    
    return target_interaction
    
    

if __name__ == '__main__':        
    cnv_file = argv[1]
    interaction_file = argv[2]
    outname = argv[3]
    if argv[4]:
        drop_chrom = argv[4]
    else:
        drop_chrom= None
    zoom(cnv_file, interaction_file, outname , drop_chrom)
