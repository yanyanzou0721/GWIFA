# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:45:26 2022

@author: yyzou

The MIT License

Copyright (c) 2022 yyzou.
"""


import click
import os
import pickle
from os.path import exists
from zoom import zoom
from fitbin import readchr, fitbin
from FS import FS



@click.command(name="GWIFA")

@click.option("pre","-p",
              default=None,
              help="if already have zoomed target region contact")

@click.option("matrix","-m",
              default=None,
              help="Contact matrix of fixed resolution")

@click.option("organize","-o",
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
              help="str or list, chromosome excluding from analysis")

@click.option("resolution","-res",
              default=100000,
              help="resolution of input matrix")

@click.option("liner_fit","-fit",
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


def GWIFA(matrix, organize, cnv_info, chrom_length_info, drop_chrom, resolution, pre=None,liner_fit=True, spacing=3, ymin=-100, ymax=100, outdir="./", outfig="test"):
    if not exists(outdir):
        os.mkdir(outdir)
    
    if pre:
        fileo = open(pre, 'rb')
        # loading data
        pac1 = pickle.load(fileo)
        # close the file
        fileo.close()
    else:
        pac1 = zoom(cnv_info, matrix, organize,drop_chrom)
    ### pac1 is a dictionary, with 2 keys ["cnv_region","target_interaction"]
    
    
    chr_len = readchr(chrom_length_info,resolution,drop_chrom)
    cumulative_interaction_intensity = fitbin(pac1["target_interaction"],organize,pac1["cnv_region"], chr_len,resolution)
    FS(cumulative_interaction_intensity,chr_len,resolution,outfig)
    
    return " done ^_^ "

if __name__ == "__main__":
    GWIFA(pre, matrix, organize, cnv_info, chrom_length_info, drop_chrom, resolution, liner_fit, spacing, ymin, ymax, outdir, outfig)