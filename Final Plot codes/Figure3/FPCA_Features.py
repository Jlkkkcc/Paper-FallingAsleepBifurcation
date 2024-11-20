#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 17:54:45 2022

@author: junheng

This code is used to run FPCA analysis for group-averaged features


"""


#%% Data preparation / Note that different MATLAB versions might save .mat into different versions, thus use either Scipy.io or mat73 based on your needs

import mat73 as mp
import scipy.io as ii
data_dict = mp.loadmat('FPCADynamics.mat')
#data_dict = ii.loadmat('FPCADynamics.mat')

ftallmat_allnorm_nonan = data_dict["ftall_mat_allnorm_strm"]
max_ck_real = data_dict["max_ck_real"]

#%% FPCA analysis

import skfda 
from skfda.exploratory.visualization import FPCAPlot
from skfda.preprocessing.dim_reduction.feature_extraction import FPCA
from skfda.representation.basis import BSpline, Fourier, Monomial

import matplotlib
matplotlib.rc('font', family='sans-serif') 

matplotlib.rc('font', serif='Helvetica Neue') 
matplotlib.rc('text', usetex='false') 
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import numpy as np
import math


# ftallmat_allnorm_nonan = ftallmat_allnorm[:,~np.isnan(ftallmat_allnorm).any(axis=0)]       #Drop column NaNs

real_len = np.size(ftallmat_allnorm_nonan,1)
epc_len = 6
epc_jmp = 3
post_epcn = math.floor(10*60/epc_len);         # Post-asleep epoch numbers

slp_onset_id = int(max_ck_real)-1

preslp = int(np.floor((60*60 -epc_len)/epc_jmp))

ftallmat_allnorm_nonan_preslp = ftallmat_allnorm_nonan[:,slp_onset_id-preslp:slp_onset_id] #Clustering only takes pre-asleep periods

grid_p = np.linspace(-(((preslp-1)*epc_jmp+epc_len)/60),0,preslp)         #Time data grids

fd = skfda.FDataGrid(
data_matrix=ftallmat_allnorm_nonan_preslp,
grid_points=grid_p,
)


#%% FPCA with L2 regularisation

from skfda.misc.regularization import L2Regularization

reg = L2Regularization()

n_cmp = 2;   #Number of components

fpca_reg = FPCA(n_components=n_cmp,regularization=reg)

pc_score_reg = fpca_reg.fit_transform(fd)  #Evaluate a score for each component time-series
fpca_reg.components_.plot()

fig11 = plt.figure(figsize=(6, 2 * 4))
FPCAPlot(
    fd.mean(),
    fpca_reg.components_,
    10,
    fig=fig11,
    n_rows=2,
).plot()

top_ft_pos = np.argmax(pc_score_reg,0)

# Mapping direclty from FPCA to feature time-series
pc_fpca_L2reg = np.squeeze(fpca_reg.components_.data_matrix) 

ft_pos_fpcaL2 = np.argmax(np.transpose(pc_score_reg),1)
print('Positively-correlated features: ', ft_pos_fpcaL2+1)
ft_neg_fpcaL2 = np.argmin(np.transpose(pc_score_reg),1)
print('Negatively-correlated features: ', ft_neg_fpcaL2+1)

#%% Figure editing

xx = fig11.get_axes()
xx[0].set_xlim([-60,0])
xx[0].set_ylim([-1.8,1.5]) 
xx[0].spines['top'].set_visible(False)
xx[0].spines['right'].set_visible(False)
xx[0].spines['bottom'].set_linewidth(3)
xx[0].spines['left'].set_linewidth(3)
xx[0].tick_params(width=3,length=6)
xx[0].set_yticks(np.linspace(-2,1.5,8)) 

xx[1].set_xlim([-60,0])
xx[1].set_ylim([-4,4]) 
xx[1].spines['top'].set_visible(False)
xx[1].spines['right'].set_visible(False)
xx[1].spines['bottom'].set_linewidth(3)
xx[1].spines['left'].set_linewidth(3)
xx[1].tick_params(width=3,length=6)

fig11






