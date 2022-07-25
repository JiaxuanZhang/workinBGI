#!/usr/bin/python
# -*- coding: utf-8 -*-
#coding=utf-8
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import sparse, stats
import h5py
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import cv2
import squidpy as sq

def main():
    """
    %prog [options]
    This softeware can be used to do saturation statistic
    """

    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--inFile", action="store", type="str", dest="inFile", help="input gene expression matrix after add img in anndata format")
    parser.add_option("-o", "--out", action="store", type="str", dest="out", help="output directory")
    parser.add_option("-f","--file", default="hello you need me", help="hello you need me")
    #parser.add_option("--tissue", action="store", type="str", dest="tissue", help="gene expression matrix only contians data under the tissue region")

    opts, args = parser.parse_args()
    if (opts.inFile == None):
        sys.exit(not parser.print_help())

    saturationFile = os.path.join(opts.out, "test_human_homology_sp.h5ad")
    os.makedirs(opts.out, exist_ok=True)
    moranl_analysis(opts.inFile, saturationFile)
    #getSaturationTable(opts.inFile, opts.tissue, saturationFile)
    #getSaturationFig(saturationFile, opts.out)


def moranl_analysis(inFile, outFile):
    #adataFile = "/jdfssz2/ST_BIOINTEL/P20Z10200N0039/06.user/liuxing2/project/clinical/mouse_EMT6/cell2loc/results/MergeRef/cell2location_map/merge_sp.h5ad"
    adata = anndata.read(inFile)
    #sc.pl.highest_expr_genes(adata, n_top=20)
    #过滤低表达基因或细胞
    sc.pp.filter_cells(adata, min_genes=3)
    #sc.pp.filter_cells(adata, max_genes=3000)
    sc.pp.filter_genes(adata, min_counts=100)
    adata.var['mt'] = adata.var_names.str.startswith(('mt-', 'MT-'))  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    #sc.pl.violin(adata, ['n_genes', 'total_counts', "pct_counts_mt"], jitter = 0.4, multi_panel=True)
    #sc.pl.scatter(adata, x='total_counts', y='n_genes')
    #归一化
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    #选择表达量前5000个的高表达基因
    highly_expr_list = adata.var['total_counts'].sort_values(ascending=False).head(5000).index.drop("Gm42418")
    adata.var.highly_expr = adata.var.index.isin(highly_expr_list)
    adata = adata[:, adata.var.highly_expr]
    #用moranl方法找高可变基因
    sq.gr.spatial_neighbors(adata)#找邻域
    from scanpy._utils import is_constant
    genes = adata.var_names.values[:4999]
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        genes=genes,
        n_perms=100,
        n_jobs=1,
    )
    genes = adata.uns["moranI"].head(2000).index.values
    adata = adata[:, adata.var.index.isin(genes)]
 #adata.write("/jdfssz2/ST_BIOINTEL/P20Z10200N0039/06.user/liuxing2/project/clinical/mouse_EMT6/cell2loc/results/MergeRef/cell2location_map/merge_human_homology_sp.h5ad")
    adata.write(outFile)

if __name__ == "__main__":
    main()
   
