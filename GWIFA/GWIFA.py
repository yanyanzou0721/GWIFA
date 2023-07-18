# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:45:26 2022

@author: yyzou

The MIT License

Copyright (c) 2022 yyzou.
"""


import click
import os
from os.path import exists
import pandas as pd
from zoom import zoom
from fitbin import readchr, fitbin
from FS import FS



@click.command(name="GWIFA")

@click.option("pre","-p",
              default=None,
              help="if already have zoomed target region contact, pass the target region .xsl")

@click.option("matrix","-m",
              default=None,
              help="Contact matrix of fixed resolution")

@click.option("organ","-o",
              default=None,
              help="Name of the subject")

@click.option("cnv_info","-cnv",
              default=None,
              help="CNV information in bed format")

@click.option("chrom_length_info","-l",
              default=None,
              help="Lenght of chromosomes, eq:hg19.len")

@click.option("drop_chrom","-d",
              default=None,
              multiple=True, 
              help="chromosomes excluding from analysis, if more than one chromosomes -d can be used for multiple times")

@click.option("resolution","-res",
              default=100000,
              help="resolution of input matrix")

@click.option("--fit/--no-fit",
              default=True,
              help="Conduct linear fit on cumulative interaction intensity")

@click.option("spacing","-s",
              default=3,
              help="spacing for second order backward difference caculation, if no linear fit is performed ")

@click.option("ymin","-ymin",
              default=-100,
              help="ylim for plot ")
@click.option("ymax","-ymax",
              default=100,
              help="ylim for plot")


@click.option("--outdir", "-O",
    default="./",
    help="path to output files.")

@click.option("--outfig", "-fig",
    default="test.pdf",
    help="name of output figure.")


def GWIFA(matrix, organ, cnv_info, chrom_length_info, drop_chrom, resolution, pre=None,fit=True, spacing=3, ymin=-100, ymax=100, outdir="./", outfig="test"):
    if not exists(outdir):
        os.mkdir(outdir)
    
    if pre:
        target_interaction = pd.read_table(pre,sep="\t")
        cnv_region = pd.read_table(cnv_info,header=None,sep="\t")
        cnv_info_txt=str(cnv_region[0].loc[0])+":"+str(cnv_region[1].loc[0])+"-"+str(cnv_region[2].loc[0])
    else:
        cnv_info_txt,target_interaction = zoom(cnv_info, matrix, outdir+organ,drop_chrom)

    tmp=[]
    for i in drop_chrom:
        tmp.append(str(i))
    drop_chrom = tmp
    print(drop_chrom)

    chr_len = readchr(chrom_length_info,resolution,drop_chrom)
    
    cumulative_interaction_intensity = fitbin(target_interaction,organ,cnv_info, chr_len,resolution)
    
    fluctuation_score,amplification_type = FS(cumulative_interaction_intensity,chr_len,resolution,outdir+outfig,ymin,ymax,fit, spacing)
    
    report_content = f"""
Genome-wide interaction fluctuation analysis
Result Report:
-------
Sample:{organ}
Amplification_region:{cnv_info_txt}

FS(fluctuation_score)={fluctuation_score}
Amplification_type={amplification_type}
    """
    with open(outdir+organ+"_"+cnv_info_txt+"_report.txt", "w") as file:
        file.write(report_content)
    
    return " done ^_^ "

if __name__ == "__main__":
    GWIFA()
