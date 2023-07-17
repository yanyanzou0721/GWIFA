# GWIFA
Genome-wide interaction fluctuation analysis for ecDNA and HSR detection

## Help info:
```
python GWIFA.py --help

```
------
## Using GWIFA from .hic file
```
1. hic2cool convert file.hic file_100k.cool -r 100000 
2. cooler dump --header --join -o file_100k.txt file_100k.cool
3. python GWIFA.py -m file_100k.txt -o fig_title -cnv cnv_region.bed -l hg19.len -d Y -d M -res 100000 -fig fig_name
```
#### Process one cnv area at a time, if there are multiple areas, set up multiple cnv_region.bed files.
------
## Input：
```
Options:
  -p TEXT              if already have zoomed target region contact( {organ}_target_region.xls ), or re-run for debug
  -m TEXT              Contact matrix of fixed resolution, full matrix generated by "cooler dump cooler dump --header --join -o output input.cool"
  -o TEXT              Name of the subject
  -cnv TEXT            CNV information in bed format
  -l TEXT              Lenght of chromosomes, eq:hg19.len
  -d TEXT              chromosomes excluding from analysis. if more than one chromosomes, -d can be used for multiple times.
  -res INTEGER         resolution of input matrix
  --fit / --no-fit     Conduct linear fit on cumulative interaction intensity
  -s INTEGER           spacing for second order backward difference
                       caculation, if no linear fit is performed
  -ymin INTEGER        ylim for plot
  -ymax INTEGER        ylim for plot
  -O, --outdir TEXT    path to output files.
  -fig, --outfig TEXT  name of output figure.
```

## Output:
    1. {organ}_target_region.xls, interaction information between  the focal amplified region and whole Genome
    2. Second_Derivation.pdf,  figure of the final result
    3. {organ}_{cnv_region}_report.txt, report of focal amplified region type and fluctuation score

## Further plot
```
Rscript circlize_plot.r {organ}_target_region.xls contact_threshold outfig
```
