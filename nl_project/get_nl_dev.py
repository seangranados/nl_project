# Purpose of this script is to quantify the deviation between experimental 
# and theoretical flux responses on the infrared detector (NIRC2) due to 
# non-linearity effects


# Some documentation:
'''
This script is designed for kp images. But can be changed for other filters.


user inputs:
        1. nirc2 epoch
        2. starfinder version which includes single psf or legacy 
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
        
    
        8. Get polynomial fit from Metchev Polynomial: y = 1.001 - 6.9*10**-6(x) - 0.70*10**-10(x**2)
        9. Get linear fit from ideal infrared detector: y = 1.001 + x
        
        10. Plot theoretical counts vs actual counts using linear fit and M Polynomial
        
        11. Get conversion of counts to % deviation formula 
        
        12. Read in PSF stars, make table to get deviation.
        
        13. Output CSV file 
        
        14. Output Plot

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

# A:
        
epoch = input('epoch (e.g. 20160503): ')
stf_ver = input('legacy or sinpsf: ')
if stf_ver == 'legacy':
    stf = 'v2_3'
if stf_ver == 'sinpsf':
    stf = 'v3_1'
    
    
# Function to identify stars not named in Grant's Work
def find_stars(ind_lis, epoch, star):
    # Steps:
    # 1. Take average of the differences between ind_frame star locations and rms_named.lis star locations
    # of known stars
    # 2. Identify this average as an offset transformation from ind_frame star locations to rms_named.lis star locations
    # 3. Apply transformation to Grant's work and match the location with the star name given in rms_named.lis
    
    # Read in rms_named.lis file, only keep named stars
    rms_named_file = f'/g/ghez/data/dr/dr1/starlists/combo/{epoch}nirc2/starfinder_{stf}/mag{epoch}nirc2_kp_rms_named.lis'
    rms_named_lis = Table.read(rms_named_file, format = 'ascii.commented_header', delimiter = '\s')
    indices_to_remove = [index for index, value in enumerate(rms_named_lis['name']) if value.startswith('star')] 
    rms_named_lis.remove_rows(indices_to_remove)
    
    # Read in ind_frame lis file, only keep named stars
    indices_to_remove = [index for index, value in enumerate(ind_lis['col1']) if value.startswith('star')]
    ind_lis.remove_rows(indices_to_remove)
    
    # Take difference between same stars in the two lists and average
    rms_named_lis.sort('name')
    ind_lis.sort('col1')
    matched_stars = []
    for s in ind_lis['col1']:
        for n in rms_named_lis['name']:
            if s == n:
                matched_stars.append(s)
            else:
                continue
    x_diff = []
    y_diff = []
    for ms in matched_stars:
        guide_rms, guide_ind = np.where(rms_named_lis['name'] == ms)[0], np.where(ind_lis['col1'] == ms)[0]
        x_diff.append(float(ind_lis['col4'][guide_ind] - rms_named_lis['x'][guide_rms]))
        y_diff.append(float(ind_lis['col5'][guide_ind] - rms_named_lis['y'][guide_rms]))
    x_rms_to_ind = np.mean(x_diff)
    y_rms_to_ind = np.mean(y_diff)
    guide = np.where(rms_named_lis['name'] == star)[0]
    x_loc = float(rms_named_lis['x'][guide]) + x_rms_to_ind
    y_loc = float(rms_named_lis['y'][guide]) + y_rms_to_ind
    return x_loc, y_loc

    
# Read in psf star list 
#psf_stf_loc = '/g/ghez/data/dr/{0}/source_list/psf_central.dat'.format('dr1')
psf_stf_loc = '/u/ciurlo/Research/starfinder_tests/PSFstars_lists/25_PSFstars/psf_central.dat'
psf_stf_lis = Table.read(psf_stf_loc, format = 'ascii.commented_header', delimiter = '\s')

# Get list of PSF stars
psf_stars = []
for psf, i in enumerate(psf_stf_lis['PSF']):
    if i == True:
        psf_stars.append(psf_stf_lis['Name'][psf])
    else:
        continue
        
# Read in psf flux star list
flux_psf_loc = '/g/ghez/data/dr/{0}/starlists/combo/{1}nirc2/starfinder_{2}/mag{1}nirc2_kp_rms_named.lis'.format('dr1', epoch, stf)

flux_psf_lis = Table.read(flux_psf_loc, format = 'ascii.commented_header', delimiter = '\s')

# Get list of mag values of psf stars
psf_mag = []
for star in psf_stars:
    guide = np.where(flux_psf_lis['name'] == star)[0]
    psf_mag.append(flux_psf_lis['m'][guide])
    
# Get zero point to convert mag to flux
zp_loc = '/g/ghez/data/dr/{0}/starlists/combo/{1}nirc2/starfinder_{2}/mag{1}nirc2_kp_0.8_stf_cal.zer'.format('dr1', epoch, stf)

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
strehl_loc = '/g/ghez/data/dr/{0}/clean/{1}nirc2/kp/strehl_source.txt'.format('dr1', epoch)
strehl_lis = Table.read(strehl_loc, format = 'ascii.no_header', delimiter = '\s')
strehl_max = np.max(strehl_lis['col2'])
guide = np.where(strehl_lis['col2'] == strehl_max)[0]
raw_frame = str((strehl_lis['col1'][guide][0])).replace('c', 'n')

raw_loc = '/g/ghez/data/dr/{0}/raw/{1}nirc2/{2}'.format('dr1', epoch, raw_frame)

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
for j, star in enumerate(psf_stars):
    for i in range(len(ind_lis['col1'])):   
        if ind_lis['col1'][i] == star:
            star_x.append(ind_lis['col4'][i])
            star_y.append(ind_lis['col5'][i])
        else:
            continue
    if (len(star_x) - 1) != j:
        star_x.append(find_stars(ind_lis, epoch, star)[0])
        star_y.append(find_stars(ind_lis, epoch, star)[1])
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






# Get poly and linear
x = np.linspace(0, 20000, 20000)
M_poly = 1.001 - 6.9*10**(-6)*(x) - 0.70*10**(-10)*(x**2)
M = 1.001 + x*M_poly
lin = 1.001 + x

# Get Counts to deviation formula
def m_poly(x):
    m_poly = 1.001 - 6.9*10**(-6)*(x) - 0.70*10**(-10)*(x**2)
    return (1.001 + x*m_poly)

def linear(x):
    linear = 1.001 + x
    return linear
    
deviation = []
for i in x:    
    actual = m_poly(i)
    theo = linear(i)
    deviation.append((1 - actual/theo)*100)

dev_list = []
x_axis = []
for counts in maxcounts:
    for i in range(len(x)):
        if int(x[i]) == int(counts):
            dev_list.append(deviation[i])
        else:
            continue
    x_axis.append(m_poly(counts))

    
dat_f = pd.DataFrame({'PSF Stars': psf_stars, 'Max Counts/COADD': maxcounts,
                   'Normalized Flux' : psf_flux,
                  'Deviation (%)': dev_list, 'x_axis': x_axis})

dat_f.to_csv('{0}/dev_{1}_{2}.csv'.format(os.getcwd(), epoch, stf_ver))
subprocess.run(['open', '{0}/dev_{1}_{2}.csv'.format(os.getcwd(), epoch, stf_ver)], check=True)

sns.set_style('darkgrid')
plt.figure(figsize = (8,8))
sns.histplot(data = dat_f, x = 'Deviation (%)')
plt.title('{0} (Frame {2}) in {1}'.format(epoch, stf_ver, raw_frame))
plt.savefig('{0}/dev_hist_{1}_{2}.jpg'.format(os.getcwd(), epoch, stf_ver), format = 'jpg')
plt.show()


# Plot 
sns.set_style('darkgrid')
plt.figure(figsize = (len(psf_stars)/1.7, len(psf_stars)/1.7))
plt.plot(M, x, label = 'Metchev Poly')
plt.plot(lin, x, label = 'Linear')
sns.scatterplot(data = dat_f, x = 'x_axis', y = 'Max Counts/COADD', hue = 'PSF Stars')
plt.xlabel('Actual Counts/COADD', fontsize = len(psf_stars)/1.5)
plt.ylabel('Theoretical Counts/COADD', fontsize = len(psf_stars)/1.5)
plt.title('{0} (Frame {2}) in {1}'.format(epoch, stf_ver, raw_frame), fontsize = len(psf_stars)/1.5)
plt.legend()
plt.savefig('{0}/dev_{1}_{2}.jpg'.format(os.getcwd(), epoch, stf_ver), format = 'jpg')
plt.show()