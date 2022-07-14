#!/usr/bin/python
import numpy as np
import pandas as pd
import scanpy as sc 
#import squidpy as sp
from anndata import AnnData
import anndata
import os, sys
from optparse import OptionParser
import csv


def main():
    """
    %prog [options]
    This softeware can be used to do saturation statistic
    """

    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--inFile", action="store", type="str", dest="inFile", help="input gene expression matrix contains reads count information")
    parser.add_option("-o", "--out", action="store", type="str", dest="out", help="output directory")
    parser.add_option("-f","--file", default="hello you need me", help="hello you need me")
    #parser.add_option("--tissue", action="store", type="str", dest="tissue", help="gene expression matrix only contians data under the tissue region")

    opts, args = parser.parse_args()
    if (opts.inFile == None):
        sys.exit(not parser.print_help())

    saturationFile = os.path.join(opts.out, "test_human_homology_sp.h5ad")
    os.makedirs(opts.out, exist_ok=True)
    m2h_homologene(opts.inFile, saturationFile)
    #getSaturationTable(opts.inFile, opts.tissue, saturationFile)
    #getSaturationFig(saturationFile, opts.out)


def m2h_homologene(inFile, outFile):
    #adataFile = "/jdfssz2/ST_BIOINTEL/P20Z10200N0039/06.user/liuxing2/project/clinical/mouse_EMT6/cell2loc/results/MergeRef/cell2location_map/merge_sp.h5ad"
    adata = anndata.read(inFile)
    human_mouse_homology_file = "/jdfssz2/ST_BIOINTEL/P20Z10200N0039/06.user/zhangjiaxuan1/mouse2human.txt"
    hmdf = pd.read_csv(human_mouse_homology_file, sep="\t")
    mouseGenes = hmdf['mouseGene'].values
    humanGenes = hmdf['humanGene'].values
    hdict = dict(zip(mouseGenes, humanGenes))
    adata.var.index = [hdict[x] if x in hdict.keys() else 'Nan' for x in adata.var.index]
    adata = adata[:, adata.var.index!='Nan']
    adata.var_names_make_unique()
    #adata.write("/jdfssz2/ST_BIOINTEL/P20Z10200N0039/06.user/liuxing2/project/clinical/mouse_EMT6/cell2loc/results/MergeRef/cell2location_map/merge_human_homology_sp.h5ad")
    adata.write(outFile)

if __name__ == "__main__":
    main()
   
