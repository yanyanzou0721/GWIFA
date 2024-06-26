# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:55:46 2022
updated on Sat Jan 30 21:29 2024
@author: yyzou
"""
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from scipy.interpolate import UnivariateSpline


def GWIFA_plot(FS,CII,CIIs,CII_2d,x_range,chr_len,title,res,ymin,ymax,outfig_name):
    
    ## plot
    plt.figure(figsize=(20,10))
    plt.suptitle(title, fontsize=30)
    
    ax = plt.subplot(1,2,1)
    x = CII["bin"].to_list()
    y = CII["cumsum"].to_list()
    ax.plot(x,y,'ro',label = 'CII')
    ax.plot(x_range,CIIs(x_range))     ####plot cumulative_interaction_intensity distribution
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,chr_len.shape[0]+1),rotation=45,size=8)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Cumulative Interaction Intensity(CII)')    
   
    ax = plt.subplot(1,2,2)
    plt.plot(x_range,CII_2d,c="darkred")
    ax.set_ylim(ymin,ymax)
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,chr_len.shape[0]+1),rotation=45,size=8)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Second-order Difference of CII')
    
    plt.savefig(outfig_name,bbox_inches="tight")



def diff(CII,HCC,periods,chr_len,res,ymin,ymax,outfig):  
    
    x_range = np.linspace(0,chr_len["bins"].sum()-2,10000,dtype="int")
    x = CII["bin"].to_list()
    y = CII["cumsum"].to_list()
    CIIs = UnivariateSpline(x,y)                 ### smooth
    ## caculate FS score
    CII_2d = CII["cumsum"].diff(periods=periods).diff(periods=periods).fillna(0).loc[x_range]/(periods^2)
    
    nd = np.sort([200 if i>200 else i for i in abs(CII_2d)])  ## polar distribution
    nd2 = nd[int(0.9*len(x_range)):]      ## top 10%
    print([nd2.sum,nd.sum()])
    FS = nd2.sum()/nd.sum()
    
    if int(HCC) >= 0.5*len(chr_len["chr"]):
        if FS < 0.8:
            amplification_type = "ecDNA"
            title = "FS = "+ str(FS) + ", ecDNA"
    else:
        amplification_type = "HSR"
        title = "FS = "+ str(FS) + ", HSR"
        
    outfig_name = outfig+"_Second_Derivation.pdf"
    
    GWIFA_plot(FS,CII,CIIs,CII_2d,x_range,chr_len,title, res,ymin,ymax*10,outfig_name)
    return FS, amplification_type
    
    

def spline_diff(CII, HCC,chr_len, res, ymin, ymax, outfig):   
    x_range = np.linspace(1,chr_len["bins"].sum()-1,10000,dtype="int") ### generate 10000 intervals between 
    
    ## caculate FS score
    x = CII["bin"].to_list()
    y = CII["cumsum"].to_list()
    CIIs = UnivariateSpline(x,y)                 ### smooth
    CIIs_2d = CIIs.derivative(n=2)(x_range)      ### second deviation
    
    nd = np.sort([200 if i>200 else i for i in abs(CIIs_2d)])  ## polar distribution
    nd2 = nd[int(0.9*len(x_range)):]      ## top 10%
    FS = nd2.sum()/nd.sum()
    
    if int(HCC) >= 0.5*len(chr_len["chr"]):
        if FS < 0.8:
            amplification_type = "ecDNA"
            title = "FS = "+ str(FS) + ", ecDNA"
    else:
        amplification_type = "HSR"
        title = "FS = "+ str(FS) + ", HSR"
    
    outfig_name = outfig+"_spline_Second_Derivation.pdf"
    GWIFA_plot(FS,CII,CIIs,CIIs_2d,x_range,chr_len,title,res,ymin,ymax,outfig_name)
    return FS,amplification_type
 
    
def FS(CII,HCC,chr_len,res,outfig,ymin=-100,ymax=100,liner_fit=True,periods=None):
    ## HCC:counts of high contact chromsomes
    if liner_fit:
        return spline_diff(CII, HCC,chr_len,res, ymin, ymax, outfig)
    else:
        return diff(CII,HCC, periods, chr_len, res, ymin, ymax, outfig)
