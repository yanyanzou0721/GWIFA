# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:55:46 2022

@author: yyzou
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

def GWIFA_plot(FS,CII,CIIs,CII_2d,x_range,chr_len, res,ymin,ymax,outfig_name):
    
    if FS < 0.8:
        title = "FS= "+ str(FS) + ", ecDNA"
    else:
        title = "FS= "+ str(FS) + ", HSR"
    
    ## plot
    plt.figure(figsize=(20,10))
    plt.suptitle(title, fontsize=30)
    
    ax = plt.subplot(1,2,1)
    x = CII["bin"].to_list()
    y = CII["cumsum"].to_list()
    ax.plot(x,y,'ro',label = 'CII')
    ax.plot(x_range,CIIs(x_range))     #### cumulative_interaction_intensity distribution
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,chr_len.shape[0]+1),rotation=45,size=8)
    
    ax = plt.subplot(1,2,2)
    plt.plot(x_range,CII_2d,c="darkred")
    ax.set_ylim(ymin,ymax)
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,chr_len.shape[0]+1),rotation=45,size=8)

    plt.savefig(outfig_name,bbox_inches="tight")

def diff(CII,periods,chr_len, res,ymin,ymax,outfig):  
    
    x_range = np.linspace(1,chr_len["bins"].sum()-1,10000,dtype="int")
    x = CII["bin"].to_list()
    y = CII["cumsum"].to_list()
    CIIs = UnivariateSpline(x,y)                 ### smooth
    ## caculate FS score
    CII_2d = CII["cumsum"].diff(periods=periods).diff(periods=periods).loc[x_range]/(periods^2)
    
    nd = np.sort([200 if i>200 else i for i in abs(CII_2d)])  ## polar distribution
    nd2 = nd[int(0.9*len(x_range)):]      ## top 10%
    FS = nd2.sum()/nd.sum()
    
    outfig_name = outfig+"_Second_Derivation.pdf"
    
    GWIFA_plot(FS,CII,CIIs,CII_2d,x_range,chr_len, res,ymin,ymax,outfig_name)
    
    

def spline_diff(CII, chr_len, res, ymin, ymax, outfig):   
    x_range = np.linspace(1,chr_len["bins"].sum()-1,10000,dtype="int") ### generate 10000 intervals between 
    
    ## caculate FS score
    x = CII["bin"].to_list()
    y = CII["cumsum"].to_list()
    CIIs = UnivariateSpline(x,y)                 ### smooth
    CIIs_2d = CIIs.derivative(n=2)(x_range)      ### second deviation
    
    nd = np.sort([200 if i>200 else i for i in abs(CIIs_2d)])  ## polar distribution
    nd2 = nd[int(0.9*len(x_range)):]      ## top 10%
    FS = nd2.sum()/nd.sum()
    
    outfig_name = outfig+"_spline_Second_Derivation.pdf"
    GWIFA_plot(FS,CII,CIIs,CIIs_2d,x_range,chr_len, res,ymin,ymax,outfig_name)
    
 
    
def FS(CII,chr_len,res,outfig,ymin=-100,ymax=100,liner_fit=True,periods=None):
    if liner_fit:
        spline_diff(CII, chr_len, res, ymin, ymax, outfig)
    else:
        diff(CII, periods, chr_len, res, ymin, ymax, outfig)
    return None