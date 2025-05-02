# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 21:57:34 2025

@author: Junheng Li/Tianyu Wei
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.io as sio

# scikit-fda (skfda) imports
import skfda
from skfda.representation import FDataGrid
from skfda.preprocessing.dim_reduction import FPCA
from skfda.misc.regularization import L2Regularization
# Optional: from skfda.exploratory.visualization import FPCAPlot

###############################################################################
# 1. Helper function(s)
###############################################################################

def innernorm(X):
    """
    This is a self-defined function for doing inner-product normalization.
    Input matrix: 
        - The rows should be the time-series for normalization.
    Output:
        - Xnorm: row-wise normalized array.
    """
    # X_inn[i] = sqrt of inner product of row i with itself
    X_inn = np.sqrt(np.diag(np.matmul(X, X.T)))
    nobs = np.size(X, 0)
    Xnorm = np.zeros_like(X)
    for i in range(nobs):
        Xnorm[i, :] = X[i, :] / X_inn[i]
    return Xnorm

###############################################################################
# 2. Load data
###############################################################################
# Load the .mat file
matfile_path = 'Ft50Dynamics_IndivFPCA.mat'
data_dict = sio.loadmat(matfile_path, struct_as_record=False, squeeze_me=True)

# Extract variables from the loaded dictionary
ftall_mat_allnorm_strm_all_subj = data_dict['ftall_mat_allnorm_strm_all_subj']
idx_subj_used = data_dict['idx_subj_used']
max_ck_real_all_subj = data_dict['max_ck_real_all_subj']

# -------------------------------
# Convert them to more convenient Python forms:
# -------------------------------

if ftall_mat_allnorm_strm_all_subj.ndim == 2 and ftall_mat_allnorm_strm_all_subj.shape[0] == 1:
    ftall_mat_allnorm_strm_all_subj = ftall_mat_allnorm_strm_all_subj[0]

if idx_subj_used.ndim == 2 and idx_subj_used.shape[0] == 1:
    idx_subj_used = idx_subj_used[0] 

if max_ck_real_all_subj.ndim == 2 and max_ck_real_all_subj.shape[0] == 1:
    max_ck_real_all_subj = max_ck_real_all_subj[0]

###############################################################################
# 3. Prepare containers to store results for all subjects
###############################################################################

all_pc_score_reg = []  # Store PC scores for each subject, shape: (n_timepoints, n_components)
all_var_exp = []       # Store variance explained, shape: (n_components,)
all_FPC1 = []          # Store FPC1 as a 1D array (size depends on subject data length)
all_FPC2 = []          # Store FPC2 as a 1D array

###############################################################################
# 4. Main loop over subjects
###############################################################################

def run_fpca_across_subjects(ftall_mat_allnorm_strm_all_subj, max_ck_real_all_subj):
    """
    Runs the FPCA routine for each subject, collects and returns:
      - pc_score_reg
      - var_exp
      - FPC1
      - FPC2
    """

    for i_subj in range(len(ftall_mat_allnorm_strm_all_subj)):
        print(i_subj)
        # ---------------------------------------------------------------------
        # (a) Load subject-specific data
        # ---------------------------------------------------------------------
        ftallmat_allnorm_nonan = ftall_mat_allnorm_strm_all_subj[i_subj]
        slp_onset_id = int(max_ck_real_all_subj[i_subj])

        # ---------------------------------------------------------------------
        # (b) Construct the time grid
        # ---------------------------------------------------------------------
        epc_len = 6
        epc_jmp = 3
        # post_epcn = math.floor(10*60/epc_len)  # You may or may not need this

        # Pre-sleep data slice
        ftallmat_allnorm_nonan_preslp = ftallmat_allnorm_nonan[:, 0:slp_onset_id]

        # Time data grid
        grid_p = np.linspace(
            -(((slp_onset_id - 1)*epc_jmp + epc_len)/60),
             0,
             slp_onset_id
        )

        # ---------------------------------------------------------------------
        # (c) Create FDataGrid and FPCA
        # ---------------------------------------------------------------------
        fd = FDataGrid(
            data_matrix=ftallmat_allnorm_nonan_preslp,
            grid_points=grid_p,
        )

        reg = L2Regularization()  # L2 regularization
        n_cmp = 2                 # Number of components

        fpca_reg = FPCA(
            n_components=n_cmp,
            regularization=reg
        )

        # Fit transform
        pc_score_reg = fpca_reg.fit_transform(fd)  # shape = (n_obs, n_components)
        var_exp = fpca_reg.explained_variance_ratio_  # shape = (n_components,)

        # ---------------------------------------------------------------------
        # (d) Extract functional principal components
        # ---------------------------------------------------------------------
        aa = fpca_reg.components_
        FPC1 = aa.data_matrix[0, :, 0]
        FPC2 = aa.data_matrix[1, :, 0]

        # ---------------------------------------------------------------------
        # (e) Accumulate results
        # ---------------------------------------------------------------------
        all_pc_score_reg.append(pc_score_reg)
        all_var_exp.append(var_exp)
        all_FPC1.append(FPC1)
        all_FPC2.append(FPC2)

    # End of for-loop

###############################################################################
# 5. Run the function
###############################################################################

# run_fpca_across_subjects will fill global lists all_pc_score_reg, etc.
run_fpca_across_subjects(ftall_mat_allnorm_strm_all_subj, max_ck_real_all_subj)

###############################################################################
# 6. Save all results in a .mat file so MATLAB can load them
###############################################################################


out_dict = {
    'all_pc_score_reg': all_pc_score_reg,  # List of arrays
    'all_var_exp': all_var_exp,           # List of arrays
    'all_FPC1': all_FPC1,                 # List of arrays
    'all_FPC2': all_FPC2                  # List of arrays
}

sio.savemat('FPCA_results_all_subjects.mat', out_dict)

print("All subjects' results saved to FPCA_results_all_subjects.mat")
