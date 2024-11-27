# Script for finding residual median for StarFinder -- Originally created by Abhimat Gautam
# Added Non-linear correction analysis -- Sean Granados
# GCG

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.io import fits
from astropy.table import Table
import pandas as pd
from scipy.stats import norm
import matplotlib

'''
Documentation:

- This script analysis the residuals from starfinder between 
3 versions (Legacy, Single PSF, and Single PSF + Non-Linear correction).

- It essentially follows the following steps:
1. Read in different data files to use 
2. take a sig cut to clean up results
3. plot residuals
4. get median of residuals

'''

# Directory of out data
out_dir = '/g/ghez/seanz/projects/single_vs_legacy/'

# Directory of pre-corrected data
combo_dir = '/g/ghez/data/dr/dr1'

# Directory of post-corrected data
nl_dir = '/g/ghez/data/dr/dr2'

stf_sin = ['v3_1','v3_2']
stf_leg = ['v2_3','v2_4']

# Combo files
with open('/g/ghez/seanz/projects/single_vs_legacy/pm_epochs.txt', 'rt') as file:
    epoch_list = file.read().split(',')
    
with open(out_dir + 'res_med_analysis.txt', 'wt') as res_med_analysis:
    res_med_analysis.write('#epoch' + '\t' + 'post_sin_median' + '\t' + 'post_sin_std' + '\t' + 'pre_sin_median' + '\t' + 'pre_sin_std' + '\t' + 'pre_leg_median' + '\t' + 'pre_leg_std' + '\t' + 'post_leg_median' + '\t' + 'post_leg_std' + '\n')

for e in epoch_list:

    e_test_pre = '{0}/combo/{1}nirc2/mag{1}nirc2_kp.fits'.format(combo_dir, e)
    e_test_post = '{0}/combo/{1}nirc2/mag{1}nirc2_kp.fits'.format(nl_dir, e)

    # Residual Files
    resid_sin_pre = '{0}/starlists/combo/{1}nirc2/starfinder_v3_1/mag{1}nirc2_kp_res.fits'.format(combo_dir, e)
    resid_sin_post = '{0}/starlists/combo/{1}nirc2/starfinder_v3_2/mag{1}nirc2_kp_res.fits'.format(nl_dir, e)
    resid_leg_pre = '{0}/starlists/combo/{1}nirc2/starfinder_v2_3/mag{1}nirc2_kp_res.fits'.format(combo_dir, e)
    resid_leg_post = '{0}/starlists/combo/{1}nirc2/starfinder_v2_4/mag{1}nirc2_kp_res.fits'.format(nl_dir, e)
    # Sig Files
    sig_pre = '{0}/combo/{1}nirc2/mag{1}nirc2_kp_sig.fits'.format(combo_dir, e)
    sig_post = '{0}/combo/{1}nirc2/mag{1}nirc2_kp_sig.fits'.format(nl_dir, e)
    # data_combo combo 
    # data_post residual post
    # data_pre residual pre
    # data_leg residual legacy

    # read in files

    #COMBO
    combo_pre = fits.open(e_test_pre)
    data_combo_pre = combo_pre[0].data
    
    combo_post = fits.open(e_test_post)
    data_combo_post = combo_post[0].data
    
    #SIG
    sig_pr = fits.open(sig_pre)
    data_sig_pr = sig_pr[0].data
    
    sig_po = fits.open(sig_post)
    data_sig_po = sig_po[0].data
    #STARFINDER
    r_sin_post = fits.open(resid_sin_post)
    data_sin_post = r_sin_post[0].data

    r_sin_pre = fits.open(resid_sin_pre)
    data_sin_pre = r_sin_pre[0].data

    r_leg_pre = fits.open(resid_leg_pre)
    data_leg_pre = r_leg_pre[0].data

    r_leg_post = fits.open(resid_leg_post)
    data_leg_post = r_leg_post[0].data

    # make 0.97 cut to remove outlier data points
    maxsigvalue_pr = np.max(data_sig_pr)
    maxsigvalue_po = np.max(data_sig_po)

    sig_cut_pr = np.where(data_sig_pr >= 0.97 * maxsigvalue_pr)
    sig_cut_po = np.where(data_sig_po >= 0.97 * maxsigvalue_po)



    # apply cut to all data arrays
    data_combo_pre = data_combo_pre[sig_cut_pr]
    data_combo_post = data_combo_post[sig_cut_po]
    data_sin_post = data_sin_post[sig_cut_po]
    data_sin_pre = data_sin_pre[sig_cut_pr]
    data_leg_post = data_leg_post[sig_cut_po]
    data_leg_pre = data_leg_pre[sig_cut_pr]

    # total pixels

    total_pixels = len(data_combo_post.flatten())


    post_sin_median = np.nanmedian(np.sqrt(data_sin_post.flatten()**2))
    post_sin_std = np.nanstd(np.sqrt(data_sin_post.flatten()**2))/np.sqrt(len(data_sin_post.flatten()))
    
    pre_sin_median = np.nanmedian(np.sqrt(data_sin_pre.flatten()**2))
    pre_sin_std = np.nanstd(np.sqrt(data_sin_pre.flatten()**2))/np.sqrt(len(data_sin_pre.flatten()))
    
    pre_leg_median = np.nanmedian(np.sqrt(data_leg_pre.flatten()**2))
    pre_leg_std = np.nanstd(np.sqrt(data_leg_pre.flatten()**2))/np.sqrt(len(data_leg_pre.flatten()))
    
    post_leg_median = np.nanmedian(np.sqrt(data_leg_post.flatten()**2))
    post_leg_std = np.nanstd(np.sqrt(data_leg_post.flatten()**2))/np.sqrt(len(data_leg_post.flatten()))
    

    
    with open(out_dir + 'res_med_analysis.txt', 'at') as res_med_analysis:
        res_med_analysis.write(e + '\t' + str(post_sin_median) + '\t' + str(post_sin_std) + '\t' + str(pre_sin_median) + '\t' + str(pre_sin_std) + '\t' + str(pre_leg_median) + '\t' + str(pre_leg_std) + '\t' + str(post_leg_median) + '\t' + str(post_leg_std) + '\n')
        
