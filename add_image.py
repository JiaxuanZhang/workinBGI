#!/usr/bin/python
# -*- coding: utf-8 -*-
#coding=utf-8
import os, sys
os.environ['OPENCV_IO_MAX_IMAGE_PIXELS']="1000000000000"
os.environ['OPENBLAS_NUM_THREADS'] = '1'
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
import tifffile as tif
import squidpy as sq


def main():
    """
    %prog [options]
    This softeware can be used to do saturation statistic
    """

    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--inFile", action="store", type="str", dest="inFile", help="input gene expression matrix in anndata format")
    parser.add_option("-img", "--inImg", action="store", type="str", dest="inImg", help="input SSDNA_img in tiff")
    parser.add_option("-o", "--out", action="store", type="str", dest="out", help="output directory")
    parser.add_option("-f","--file", default="hello you need me", help="hello you need me")
    #parser.add_option("--tissue", action="store", type="str", dest="tissue", help="gene expression matrix only contians data under the tissue region")

    opts, args = parser.parse_args()
    if (opts.inFile == None):
        sys.exit(not parser.print_help())

    saturationFile = os.path.join(opts.out, "add_image_sp.h5ad")
    os.makedirs(opts.out, exist_ok=True)
    add_img(opts.inFile,opts.inImg, saturationFile)
    #getSaturationTable(opts.inFile, opts.tissue, saturationFile)
    #getSaturationFig(saturationFile, opts.out)


def add_img(inFile,inImg, outFile):
    #adataFile = "/jdfssz2/ST_BIOINTEL/P20Z10200N0039/06.user/liuxing2/project/clinical/mouse_EMT6/cell2loc/results/MergeRef/cell2location_map/merge_sp.h5ad"
    adata = anndata.read(inFile)
    img = tif.imread(inImg)
    adata.uns['spatial'] = dict()
    adata.uns['spatial']['library_id'] = dict() #library_id可改为芯片名
    adata.uns['spatial']['library_id']['images'] = dict()
    adata.uns['spatial']['library_id']['scalefactors']= dict()
    adata.uns['spatial']['library_id']['scalefactors']['tissue_lowres_scalef']= 0.051635113
    adata.uns['spatial']['library_id']['scalefactors']['tissue_hires_scalef']= 0.17211704
    adata.uns['spatial']['library_id']['scalefactors']['spot_diameter_fullres']= 89.50991269308707
    adata.uns['spatial']['library_id']['scalefactors']['fiducial_diameter_fullres']= 144.59293588883295
    img2 = img[::6, ::6] #图像缩小比例，需自调，暂设为6
    adata.uns['spatial']['SS200000681TL_A1']['images']['hires'] = img2
    #adata.obsm['spatial']=adata.obsm['spatial']+20
    #adata.obsm['spatial'][:,1] = adata.obsm['spatial'][:,1]+200 #上层点移动的位置以适应底图，需自调
    #sc.pl.spatial(adata[~adata.obs['leiden'].isin(['2']), :], img_key="hires", color=['leiden'], save=[Ture],spot_size=50) #去除边缘的2号聚类点，保存图像，默认路径保存
    #adata.write("/jdfssz2/ST_BIOINTEL/P20Z10200N0039/06.user/liuxing2/project/clinical/mouse_EMT6/cell2loc/results/MergeRef/cell2location_map/merge_human_homology_sp.h5ad")
    adata.write(outFile)

if __name__ == "__main__":
    main()
   
