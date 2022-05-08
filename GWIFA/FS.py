# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:55:46 2022

@author: yyzou
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from scipy.interpolate import UnivariateSpline

def diff(CII,periods,chr_len, res,ymin,ymax,outfig):  
    
    x_range = np.linspace(1,chr_len["bins"].sum(),10000,dtype="int")
    
    plt.figure(figsize=(10,10))
    plt.suptitle(outfig, fontsize=30)
    """
    ax = plt.subplot(1,2,1)
    x = CII["case"]["bin"].to_list()
    y = CII["case"]["cumsum"].to_list()
    CIIs_case = UnivariateSpline(x,y)
    #ax.semilogy(x,y,'ro',label = 'case')
    #ax.semilogy(x_range,CIIs_case(x_range))
    ax.plot(x,y,'ro',label = 'case')
    ax.plot(x_range,CIIs_case(x_range))
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,24),rotation=45,size=8)
    """         
    ax = plt.subplot(1,1,1)
    plt.plot(x_range,CII["case"]["cumsum"].diff(periods=periods).diff(periods=periods).loc[x_range]/(periods^2),c="darkred")
    #plt.plot(x_range,CIIs_case.derivative(n=2)(x_range))
    ax.set_ylim(ymin,ymax)
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,24),rotation=45,size=8)
    """         
    ax = plt.subplot(2,2,2)         
    x = CII["control"]["bin"].to_list()
    y = CII["control"]["cumsum"].to_list()
    CIIs_control = UnivariateSpline(x,y)  ### 1-D smoothing spline fit to a given set of CII points.
    #ax.semilogy(x,y,'ro',label = 'control') 
    #ax.semilogy(x_range,CIIs_control(x_range))
    ax.plot(x,y,'ro',label = 'case')
    ax.plot(x_range,CIIs_case(x_range))
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,24),rotation=45,size=8)
    
    ax = plt.subplot(2,2,4)
    plt.plot(x_range,CII["control"]["cumsum"].diff(periods=periods).diff(periods=periods).loc[x_range]/(periods^2),c="darkred")
    #plt.plot(x_range,CIIs_control.derivative(n=2)(x_range))
    ax.set_ylim(ymin,ymax)
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,24),rotation=45,size=8)
    """
    plt.savefig(outfig+"_Second_Differences.pdf",bbox_inches="tight")

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
    
    if FS < 0.8:
        title = "FS= "+ str(FS) + ", ecDNA"
    else:
        title = "FS= "+ str(FS) + ", HSR"
        
    ## plot
    plt.figure(figsize=(20,10))
    plt.suptitle(title, fontsize=30)
    ax = plt.subplot(1,2,1)
    ax.plot(x,y,'ro',label = 'CII')
    ax.plot(x_range,CIIs(x_range))     #### cumulative_interaction_intensity distribution
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,chr_len.shape[0]+1),rotation=45,size=8)
             
    ax = plt.subplot(1,2,2)
    plt.plot(x_range,CIIs_2d,c="darkred")
    ax.set_ylim(ymin,ymax)
    ax.set_xticks(chr_len["bin_num"].to_list())
    ax.set_xticklabels(range(1,chr_len.shape[0]+1),rotation=45,size=8)

    plt.savefig(outfig+"_spline_Second_Derivation.pdf",bbox_inches="tight")
    
 
    
def FS(CII,chr_len,res,outfig,ymin=-100,ymax=100,liner_fit=True,periods=None):
    if liner_fit:
        spline_diff(CII, chr_len, res, ymin, ymax, outfig)
    else:
        diff(CII, periods, chr_len, res, ymin, ymax, outfig)
    return None