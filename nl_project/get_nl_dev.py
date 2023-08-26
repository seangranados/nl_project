# Purpose of this script is to quantify the deviation between experimental 
# and theoretical flux responses on the infrared detector (NIRC2) due to 
# non-linearity effects


# Some documentation:
'''
This script is designed for kp images. But can be changed for other filters.

user inputs:
        1. data release version: dr1 or dr2 (non-linearity correction + dark subtraction)
        2. nirc2 epoch that exists in dr1 or dr2
        3. starfinder version which includes single psf or legacy 
steps:
        1. obtain list of string names for psf stars used in starfinder 
        from the psf_central.dat file located in 
        /g/ghez/data/dr/{data release version}/source_list/psf_central.dat
        
        2. obtain magnitudes of psf stars located in
        /g/ghez/data/dr/{drversion}/starlists/combo/{epoch}nirc2/starfinder_{stfversion}/mag{epoch}nirc2_kp_rms_named.lis
        
        3. obtain zero point approximation located in 
        '/g/ghez/data/dr/{0}/starlists/combo/{1}nirc2/starfinder_{2}/mag{1}nirc2_kp_0.8_stf_cal.zer'

        4. find normalized flux for psf stars by using standard mag to flux formula
        (normalized to S3-22)
        
        5. find max counts / coadd from highest strehl raw frame for a given epoch. This is 
        calculated from '/g/ghez/data/dr/{0}/clean/{1}nirc2/kp/strehl_source.txt'
        
        6. Use Grant's individual frame stf runs to locate psf stars in raw frame. This is 
        located in 
        '/g/ghez/data/dr/dr1/starlists/ind_frames/{0}nirc2/kp/starfinder_v3_0_0/align*/lis/{1}_0.6_stf_cal.l
        
        7. locate psf star in raw frame, cut out 80x80 pixel box (arbitrary), find max count value, 
        divide by 10 (for kp)
        
        8. with normalized flux and maxcounts/coadd, plot
        
        9. find deviation from a linear fit (this is done through the first two stars, 
        some epochs have the second star lower than the first, resulting in an incorrect linear
        trend)
        
        10. Output figure and csv file with relevant values

'''



# Import modules
import numpy as np
import pandas as pd
import heapq
import os
import glob
import seaborn as sns
import subprocess
import warnings
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

warnings.filterwarnings("ignore")

mode = input('dr1 or dr2: ')
epoch = input('epoch (e.g. 20160503): ')
stf_ver = input('legacy or sinpsf: ')
if stf_ver == 'legacy' and mode == 'dr1':
    stf = 'v2_3'
if stf_ver == 'legacy' and mode == 'dr2':
    stf = 'v2_4'
if stf_ver == 'sinpsf' and mode == 'dr1':
    stf = 'v3_1'
if stf_ver == 'sinpsf' and mode == 'dr2':
    stf = 'v3_2'
    
# Read in psf star list 
psf_stf_loc = '/g/ghez/data/dr/{0}/source_list/psf_central.dat'.format(mode)

psf_stf_lis = Table.read(psf_stf_loc, format = 'ascii.commented_header', delimiter = '\s')

# Get list of PSF stars
psf_stars = []
for psf, i in enumerate(psf_stf_lis['PSF']):
    if i == True:
        psf_stars.append(psf_stf_lis['Name'][psf])
    else:
        continue
        
# Read in psf flux star list
flux_psf_loc = '/g/ghez/data/dr/{0}/starlists/combo/{1}nirc2/starfinder_{2}/mag{1}nirc2_kp_rms_named.lis'.format(mode, epoch, stf)

flux_psf_lis = Table.read(flux_psf_loc, format = 'ascii.commented_header', delimiter = '\s')

# Get list of mag values of psf stars
psf_mag = []
for star in psf_stars:
    guide = np.where(flux_psf_lis['name'] == star)[0]
    psf_mag.append(flux_psf_lis['m'][guide])
    
# Get zero point to convert mag to flux
zp_loc = '/g/ghez/data/dr/{0}/starlists/combo/{1}nirc2/starfinder_{2}/mag{1}nirc2_kp_0.8_stf_cal.zer'.format(mode, epoch, stf)

with open(zp_loc, 'rt') as table_zp_loc:
    zp = table_zp_loc.readlines()

# Get zero point
zp = float(str(zp[4])[5:10])

# Get list of normalized flux values to S3-22
psf_flux = []
guide = np.where(flux_psf_lis['name'] == 'S3-22')[0]
s3_22_mag = flux_psf_lis['m'][guide]
s3_22 = 10**((zp-s3_22_mag)/2.5)
for mag in psf_mag:
    psf_flux.append(float((10**((zp-mag)/2.5))/s3_22))
    
# Get max counts/coadd of highest strehl raw frame of image
strehl_loc = '/g/ghez/data/dr/{0}/clean/{1}nirc2/kp/strehl_source.txt'.format(mode, epoch)
strehl_lis = Table.read(strehl_loc, format = 'ascii.no_header', delimiter = '\s')
strehl_max = np.max(strehl_lis['col2'])
guide = np.where(strehl_lis['col2'] == strehl_max)[0]
raw_frame = str((strehl_lis['col1'][guide][0])).replace('c', 'n')

raw_loc = '/g/ghez/data/dr/{0}/raw/{1}nirc2/{2}'.format(mode, epoch, raw_frame)

hdu = fits.open(raw_loc, ignore_missing_end = True)
wcs = WCS(hdu[0].header)
data = hdu[0].data

star_loc = []

# Need individual frame locations, use Grant's work 
ind_loc = '/g/ghez/data/dr/dr1/starlists/ind_frames/{0}nirc2/kp/starfinder_v3_0_0/align*/lis/{1}_0.6_stf_cal.lis'.format(epoch,(str((strehl_lis['col1'][guide][0])))[0:5])
match = (glob.glob(ind_loc))[0]


ind_lis = Table.read(match, format = 'ascii.no_header', delimiter = '\s')
star_x = []
star_y = []
for star in psf_stars:
    for i in range(len(ind_lis['col1'])):   
        if ind_lis['col1'][i] == star:
            star_x.append(ind_lis['col4'][i])
            star_y.append(ind_lis['col5'][i])
        else:
            continue
            
for i in range(len(psf_stars)):
    star_loc.append((star_x[i], star_y[i]))
    
star_size = (80,80) # arbitrary, feel free to change

cutout = []
for i in range(len(star_loc)):
    cutout.append((Cutout2D(data, star_loc[i], star_size, wcs=wcs)).data)
    
maxcounts = []
for i in range(len(cutout)):
    maxcounts.append(np.max(cutout[i])/10)
    
# Get linear data by plotting linear fit going through first two stars
p1 = ((np.min(psf_flux)), (np.min(maxcounts)))
p2 = ((heapq.nsmallest(2, psf_flux)[1]),(heapq.nsmallest(2, maxcounts)[1]))

x = [p1[0],p2[0]]
y = [p1[1],p2[1]]

slope, intercept = np.polyfit(x, y, 1)

x_line = np.linspace(min(x), max(psf_flux), 100)
y_line = slope * x_line + intercept

# Get deviation
dev = []
for i, mc in enumerate(maxcounts):
    dev.append(np.abs(1 - (mc/(slope*psf_flux[i] + intercept)))*100)

# Output csv file using dataframe and plot
df = pd.DataFrame({'PSF Stars': psf_stars, 'Max Counts/COADD': maxcounts,
                   'Normalized Flux' : psf_flux,
                  'Deviation (%)': dev})

df.to_csv('{0}/dev_{1}_{2}.csv'.format(os.getcwd(), epoch, mode))
subprocess.run(['open', '{0}/dev_{1}_{2}.csv'.format(os.getcwd(), epoch, mode)], check=True)

sns.set_style('darkgrid')
plt.figure()
sns.scatterplot(data = df, x = 'Normalized Flux', y = 'Max Counts/COADD', hue = 'PSF Stars')
plt.plot(x_line, y_line)
plt.xlim(0,7)
plt.ylim(0,20000)
plt.title('{0} in {1} for Frame {2}: Deviation'.format(epoch, mode, raw_frame))
plt.savefig(os.getcwd() + '/dev_{0}_{1}.jpg'.format(epoch, mode), format = 'jpg')
plt.show()