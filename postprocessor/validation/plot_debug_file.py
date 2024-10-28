#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 15:32:51 2021

@author: paolo
"""
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--debug_file', required=True,
                    help="Input debug file")
parser.add_argument('-o', '--outdir', default=None,
                    help="Output directory")
parser.add_argument('-s', '--plot_selected_channels', type=bool, default=False,
                    help="Boolean flag for plotting selected channels")
args = parser.parse_args()

data = Dataset(args.debug_file,'r')

ref_x = data['Sounder_FOV_Wavenumbers'][:]
ref_y = data['Sounder_FOV_Radiance'][:]
fm_x = data['FM_Wavenumbers'][:]
fm_y  = data['FM_Radiance_1'][:]

figsize=(18,13)

if args.plot_selected_channels:
    ind = data['FM_Indices'][:]
    fm_x = fm_x[ind]
    ref_x = ref_x[ind]
    ref_y = ref_y[ind]
    fm_y  = fm_y[ind]
plt.figure(figsize=figsize)


plt.plot(ref_x,ref_y,lw=1,c='grey')
plt.plot(fm_x,fm_y,c='gold',lw=.8,alpha=0.8)


plt.scatter(ref_x,ref_y,marker='o',s=5, label='SounderFOV')
plt.scatter(fm_x,fm_y,marker='s',s=5, label='FM')
plt.xlabel(r'Wavenumbers [$cm^{-1}$]',size=14)
plt.legend(prop={'size':15})

if args.outdir is None:
    plt.show()
else:
    plt.savefig(args.outdir + '/radiances.png',dpi=300,bbox_inches='tight'      )

plt.figure(figsize=figsize)
plt.bar(ref_x,ref_y - fm_y,color='darkgreen',align='center')
plt.axhline(0,c='grey',ls='dashed',lw=1)
plt.xlabel(r'Wavenumbers [$cm^{-1}$]',size=15)
plt.ylabel(r'$\Delta Rad$',size=15)

if args.outdir is None:
    plt.show()
else:
    plt.savefig(args.outdir + '/radiance_diff.png',dpi=300,bbox_inches='tight')

