#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import subprocess
import shutil
import numpy as np
from astropy.table import Table
import pandas as pd
from astropy.table import vstack
from astropy.io import fits

class Lp_analyze(object):

    # Constructor
    def __init__(self, epoch_lp, epoch_kp, corr):

        # Attributes
        
        # Epochs to work on
        self.epoch_lp = epoch_lp
        self.epoch_kp = epoch_kp

        # Location of SgrA in pixels of Lp epoch
        self.sgra_file = pd.read_csv('/g/ghez/seanz/projects/G-objects/orbital_parameters/sgra_positions_v4_1_04.csv')

        # Extract X and Y directly in a single line
        self.sgra_loc = self.sgra_file.loc[self.sgra_file['Epoch'].astype(str) == self.epoch_lp[:8], ['X [pixels]', 'Y [pixels]']].iloc[0]


        # Region of interest
        self.xpixels = 100
        self.ypixels = 100
        
        
        # Default locations for data, assuming starfinder has been run with newest version
        self.lp_dir = '/g/ghez/seanz/projects/G-objects/data/NIRC2/'
        self.kp_dir = '/g/ghez/data/dr/dr2/'
        self.orbits_dir = '/g/ghez/seanz/source_list/orbits.dat'
        self.label_dir = '/g/ghez/seanz/source_list/label.dat'

        # Align flags 
        self.alignFlags = '-R 3 -p -a 0'

        # StarFinder version
            # 4_1: newest version of StarFinder for Lp (e.g. fixpos bug fix, non-linearity corrected, 90px PSF FWHM size, 0.7 corrMain)
            # where the calibrator stars differ from the usual GCG Kp choices (S0-12, S1-1, S1-20, S1-17)

            # 4_0: newest version for StarFinder for Kp (e.g. fixpos bug fix, non-linearity corrected, 170px PSF FWHM size, 0.8 corrMain)
            # where the calibrator stars are the usual GCG Kp choice (irs16NW, S3-22, S1-17, S1-34, S4-3, S1-1, S1-21, S3-370, S3-88, S3-36, S2-63)
        
        self.stf_lp_ver = '4_1_04'
        self.stf_kp_ver = '4_1'


        self.stf_lp_lis = ['0.4_stf', '1_0.4_stf', '2_0.4_stf', '3_0.4_stf']
        self.stf_kp_lis = ['0.8_stf', '1_0.6_stf', '2_0.6_stf', '3_0.6_stf']

        self.working_dir = '/g/ghez/seanz/projects/G-objects/starplanter/'


        # Flags for force mode
        self.deblend = "0"
        self.makePsf = "0"
        self.makeRes = "0"
        self.makeStars = "0"
        self.cooStar = "irs16C"
        self.corr = corr
        self.starlist = "/g/ghez/seanz/source_list/psf_central.dat"
        self.trimfake = "0"
        self.force = "1"
        self.fix_pos = "1"
        self.legacy = "0"

        
        self.calibrators = ["S3-22", "S1-17", "S4-6", "S4-3", "S3-370", "S3-88"]
        # self.calibrators = ["S1-20", "S0-12", "S1-1", "S1-17"]




    def sgra(self):
        return self.sgra_loc



    def setup_align(self):
        """Make align directory, write align.list files, and copy over starlists"""
        
        # Make a directory to store align outputs
        os.makedirs(f'align/{self.epoch_lp}_v{self.stf_lp_ver}', exist_ok = True)
        shutil.copy(self.orbits_dir, f'align/{self.epoch_lp}_v{self.stf_lp_ver}/')
        shutil.copy(self.label_dir, f'align/{self.epoch_lp}_v{self.stf_lp_ver}/')

        
        # Write Main map align.list file
        with open(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[0][0:3]}_named.list', 'w') as align:
            align.write(f'mag{self.epoch_lp}_lp_{self.stf_lp_lis[0]}_cal.lis 8 ref' + '\n')
            align.write(f'mag{self.epoch_kp}_kp_{self.stf_kp_lis[0]}_cal.lis 8')
    
        # Copy over cal.lis files
        shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/mag{self.epoch_lp}_lp_{self.stf_lp_lis[0]}_cal.lis', f'align/{self.epoch_lp}_v{self.stf_lp_ver}/')
        shutil.copy(self.kp_dir + f'starlists/combo/{self.epoch_kp}/starfinder_v{self.stf_kp_ver}/mag{self.epoch_kp}_kp_{self.stf_kp_lis[0]}_cal.lis', f'align/{self.epoch_lp}_v{self.stf_lp_ver}/')


        # Sub maps
        for i in range(1,3+1):
            with open(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[i][0:5]}_named.list', 'w') as align:
                align.write(f'm{self.epoch_lp}_lp_{self.stf_lp_lis[i]}_cal.lis 8 ref' + '\n')
                align.write(f'm{self.epoch_kp}_kp_{self.stf_kp_lis[i]}_cal.lis 8')
    
            # Copy over cal.lis files
            shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/m{self.epoch_lp}_lp_{self.stf_lp_lis[i]}_cal.lis', f'align/{self.epoch_lp}_v{self.stf_lp_ver}/')
            shutil.copy(self.kp_dir + f'starlists/combo/{self.epoch_kp}/starfinder_v{self.stf_kp_ver}/m{self.epoch_kp}_kp_{self.stf_kp_lis[i]}_cal.lis', f'align/{self.epoch_lp}_v{self.stf_lp_ver}/')
        
        print("Aligned directories made...")

    def run_align(self):
        """Run align on main map and sub maps"""

        # Run in align directory 
        os.chdir(f'align/{self.epoch_lp}_v{self.stf_lp_ver}')
        
        # Run main map align
        cmd = f'java align {self.alignFlags} -accel_file /g/ghez/seanz/projects/G-objects/starplanter/align/{self.epoch_lp}_v{self.stf_lp_ver}/label.dat -o /g/ghez/seanz/projects/G-objects/starplanter/align/{self.epoch_lp}_v{self.stf_lp_ver}/orbits.dat -magMix -r align_kp_{self.stf_lp_lis[0][0:3]}_named align_kp_{self.stf_lp_lis[0][0:3]}_named.list'

        subprocess.call(cmd, shell=True)

        # Run sub map align
        for i in range(1,3+1):
            cmd = f'java align {self.alignFlags} -accel_file /g/ghez/seanz/projects/G-objects/starplanter/align/{self.epoch_lp}_v{self.stf_lp_ver}/label.dat -o /g/ghez/seanz/projects/G-objects/starplanter/align/{self.epoch_lp}_v{self.stf_lp_ver}/orbits.dat -magMix -r align_kp_{self.stf_lp_lis[i][0:5]}_named align_kp_{self.stf_lp_lis[i][0:5]}_named.list'

            subprocess.call(cmd, shell = True)






    def setup_force(self):
        "Make force mode directory, copy over fits files for force mode to work on"

        # Change back to starplanter directory
        os.chdir(self.working_dir)

        # Make force mode directory
        os.makedirs(f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/', exist_ok = True)


        # Main map setup

        # Copy over image to work on
        shutil.copy(f"{self.lp_dir}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp.fits", f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp.fits')


        # input -- make it the missed kp stars (align.miss0)
        shutil.copy(f"align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[0][0:3]}_named.miss0", f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp_inp.lis')


        # Copy over necessary files 
        
        # back.fits
        shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/mag{self.epoch_lp}_lp_back.fits', f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp_back.fits')
        
        # psf.fits
        shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/mag{self.epoch_lp}_lp_psf.fits', f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp_psf.fits')
        
        # .coo
        shutil.copy(self.lp_dir + f'/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp.coo', f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp.coo')
        
        




        # Sub map setup
        for i in range(1,3+1):

            # Copy over image to work on
            shutil.copy(f"{self.lp_dir}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}.fits", f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}.fits')
            
            # input -- make it the missed kp stars (align.miss0)
            shutil.copy(f"align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[i][0:5]}_named.miss0", f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}_inp.lis')

            # Copy over necessary files 
            
            # back.fits
            shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/m{self.epoch_lp}_lp_{i}_back.fits', f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}_back.fits')
            
            # psf.fits
            shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/m{self.epoch_lp}_lp_{i}_psf.fits', f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}_psf.fits')
            
            # .coo
            shutil.copy(self.lp_dir + f'/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}.coo', f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}.coo')


        print("""Forced mode directories made...
        
        """)

    def missed_kp(self):

        print("Locating missed K' sources in L' image...")
        
        """Run force mode on Lp image with the missed aligned Kp stars as a double check"""
        # Run in force mode directory
        os.chdir(f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/')
        
        # Main map run
        
        # Image for force mode to work on
        imageFile = f'{os.getcwd()}/mag{self.epoch_lp}_lp.fits'

        stf_inp_lis = f"{os.getcwd()}/mag{self.epoch_lp}_lp_inp.lis"

        fileIDLbatch = f'idlbatch_{self.epoch_lp}_' + 'lp_main'
        fileIDLlog = fileIDLbatch + '.log'
        _batch = open(fileIDLbatch, 'w', )
        _batch.write("find_stf_lis, ")
        _batch.write(f"\"{imageFile}\"" + ", ")
        _batch.write(f"{self.corr}"  + ", ")
        _batch.write(f"\"{stf_inp_lis}\"" + ", ")
        _batch.write(f"deblend={self.deblend}" + ", ")
        _batch.write(f"makePsf={self.makePsf}" + ", ")
        _batch.write(f"makeRes={self.makeRes}" + ", ")
        _batch.write(f"makeStars={self.makeStars}" + ", ")
        _batch.write(f"legacy={self.legacy}" + ", ")
        _batch.write(f"cooStar=\"{self.cooStar}\"" + ", ")
        _batch.write(f"starlist=\"{self.starlist}\"" + ", ")
        _batch.write(f"trimfake={self.trimfake}" + ", ")
        _batch.write(f"force={self.force}" + ", ")
        _batch.write(f"fix_pos={self.fix_pos}")
        _batch.write("\n")
        _batch.write("exit\n")
        _batch.close()
    
        # run forcemode from batch file
        cmd = 'idl < ' + fileIDLbatch + ' >& ' + fileIDLlog   
        print(cmd)
        subprocess.call(cmd, shell=True)


        for i in range(1,3+1):
           # Image for force mode to work on
        
            imageFile = f"{os.getcwd()}/m{self.epoch_lp}_lp_{i}.fits"
            
            stf_inp_lis = f"{os.getcwd()}/m{self.epoch_lp}_lp_{i}_inp.lis"

            fileIDLbatch = f'idlbatch_{self.epoch_lp}_' + f'lp_{i}'
            fileIDLlog = fileIDLbatch + '.log'
            _batch = open(fileIDLbatch, 'w', )
            _batch.write("find_stf_lis, ")
            _batch.write(f"\"{imageFile}\"" + ", ")
            _batch.write(f"{self.corr}"  + ", ")
            _batch.write(f"\"{stf_inp_lis}\"" + ", ")
            _batch.write(f"deblend={self.deblend}" + ", ")
            _batch.write(f"makePsf={self.makePsf}" + ", ")
            _batch.write(f"makeRes={self.makeRes}" + ", ")
            _batch.write(f"makeStars={self.makeStars}" + ", ")
            _batch.write(f"legacy={self.legacy}" + ", ")
            _batch.write(f"cooStar=\"{self.cooStar}\"" + ", ")
            _batch.write(f"starlist=\"{self.starlist}\"" + ", ")
            _batch.write(f"trimfake={self.trimfake}" + ", ")
            _batch.write(f"force={self.force}" + ", ")
            _batch.write(f"fix_pos={self.fix_pos}")
            _batch.write("\n")
            _batch.write("exit\n")
            _batch.close()
            
            # run forcemode from batch file
            cmd = 'idl < ' + fileIDLbatch + ' >& ' + fileIDLlog   
            print(cmd)
            subprocess.call(cmd, shell=True)




    def align_force_compare(self):
        """Compare the align.pos stars to the force mode stars
        we want to avoid double starplanting stars"""

        os.chdir(self.working_dir)
        
        # Align.pos files

        # read in align.pos main 
        align_pos_main = Table.read(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[0][0:3]}_named.pos', format = 'ascii.no_header', delimiter = '\s')

        
        ax_main = np.array(list(align_pos_main['col2']))
        ay_main = np.array(list(align_pos_main['col3']))

        # read in align.pos submap 1
        align_pos_1 = Table.read(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[1][0:5]}_named.pos', format = 'ascii.no_header', delimiter = '\s')

        ax_1 = np.array(list(align_pos_1['col2']))
        ay_1 = np.array(list(align_pos_1['col3']))


        # read in align.pos submap 2
        align_pos_2 = Table.read(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[2][0:5]}_named.pos', format = 'ascii.no_header', delimiter = '\s')

        ax_2 = np.array(list(align_pos_2['col2']))
        ay_2 = np.array(list(align_pos_2['col3']))


        # read in align.pos submap 3
        align_pos_3 = Table.read(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[3][0:5]}_named.pos', format = 'ascii.no_header', delimiter = '\s')

        ax_3 = np.array(list(align_pos_3['col2']))
        ay_3 = np.array(list(align_pos_3['col3']))

        align_x = [ax_main, ax_1, ax_2, ax_3]
        align_y = [ay_main, ay_1, ay_2, ay_3]

        # Force mode starlists

        # read in force mode main
        force_main = Table.read(f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp_{self.corr}_stf.lis', format = 'ascii.no_header', delimiter = '\s')

        fx_main = np.array(list(force_main['col4']))
        fy_main = np.array(list(force_main['col5']))

        
        # read in force mode submap 1
        force_1 = Table.read(f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_1_{self.corr}_stf.lis', format = 'ascii.no_header', delimiter = '\s')

        fx_1 = np.array(list(force_1['col4']))
        fy_1 = np.array(list(force_1['col5']))

        # read in force mode submap 2
        force_2 = Table.read(f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_2_{self.corr}_stf.lis', format = 'ascii.no_header', delimiter = '\s')

        fx_2 = np.array(list(force_2['col4']))
        fy_2 = np.array(list(force_2['col5']))

        # read in force mode submap 3
        force_3 = Table.read(f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_3_{self.corr}_stf.lis', format = 'ascii.no_header', delimiter = '\s')

        fx_3 = np.array(list(force_3['col4']))
        fy_3 = np.array(list(force_3['col5']))

        force_x = [fx_main, fx_1, fx_2, fx_3]
        force_y = [fy_main, fy_1, fy_2, fy_3]

        
        close_pairs = []
        for x1, y1, x2, y2 in zip(align_x, align_y, force_x, force_y):
            close_pairs.append(self.stars_match(x1, y1, x2, y2, 2))

        print(f"""---COMPARE ALIGN AND FORCE STARLISTS---

        ---MAIN---
        NUMBER OF STARS (ALIGN): {(ax_main.shape)[0]}
        NUMBER OF STARS (FORCE): {(fx_main.shape)[0]}

        MATCHED STARS ({len(close_pairs[0])}):
        ---SUB1---
        NUMBER OF STARS (ALIGN): {(ax_1.shape)[0]}
        NUMBER OF STARS (FORCE): {(fx_1.shape)[0]}

        MATCHED STARS ({len(close_pairs[1])}):
        ---SUB2---
        NUMBER OF STARS (ALIGN): {(ax_2.shape)[0]}
        NUMBER OF STARS (FORCE): {(fx_2.shape)[0]}

        MATCHED STARS ({len(close_pairs[2])}):
        ---SUB3---
        NUMBER OF STARS (ALIGN): {(ax_3.shape)[0]}
        NUMBER OF STARS (FORCE): {(fx_3.shape)[0]}

        MATCHED STARS ({len(close_pairs[3])}):
        
        
        """)



        


        print("Obtain initial starplanting list from missed kp sources...")

        # Loop over each index in reverse to remove from force_main
        for i in reversed(range(len(close_pairs[0]))):
            force_main.remove_rows(close_pairs[0][i][1])  # Access the specific index for removal
        
        # Repeat for other tables
        for i in reversed(range(len(close_pairs[1]))):
            force_1.remove_rows(close_pairs[1][i][1])
        
        for i in reversed(range(len(close_pairs[2]))):
            force_2.remove_rows(close_pairs[2][i][1])
        
        for i in reversed(range(len(close_pairs[3]))):
            force_3.remove_rows(close_pairs[3][i][1])


        
        
        print("""
        
        Save current forced stf lis files as initial starplanting list...
        
        The file structure needed for starplanting is: 
        
        name     m     x     y     flux
        
        """)


        def create_force_sp_table(force_data):
            return pd.DataFrame({
                'name': np.array(force_data['col1']),
                'm': np.array(force_data['col2']),
                'x': np.array(force_data['col4']),  # x is column 4
                'y': np.array(force_data['col5']),  # y is column 5
                'flux': np.array(force_data['col9'])
            })
        
        force_sp_main = create_force_sp_table(force_main)
        force_sp_1 = create_force_sp_table(force_1)
        force_sp_2 = create_force_sp_table(force_2)
        force_sp_3 = create_force_sp_table(force_3)

        os.makedirs(f"lp_model/{self.epoch_lp}/", exist_ok = True)
        force_sp_main.to_csv(f"lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[0]}_model.lis", sep = '\t', header = False, index = False)
        force_sp_1.to_csv(f"lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[1]}_model.lis", sep = '\t', header = False, index = False)
        force_sp_2.to_csv(f"lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[2]}_model.lis", sep = '\t', header = False, index = False)
        force_sp_3.to_csv(f"lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[3]}_model.lis", sep = '\t', header = False, index = False)
        



    def stars_match(self, x1, y1, x2, y2, threshold):
        """Matching stars between lists"""

        dx = x1[:, np.newaxis] - x2[np.newaxis, :]
        dy = y1[:, np.newaxis] - y2[np.newaxis, :]
        dist_sq = dx**2 + dy**2

        close_pairs = np.argwhere(dist_sq <= threshold**2)

        close_pairs_list = [tuple(pair) for pair in close_pairs]

        return close_pairs_list


    def kp_and_lp_sources(self):

        print("""
        
        Identifying sources in L' and K'...
        
        """)

        # Main map:
        align_pos_main = Table.read(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[0][0:3]}_named.pos', 
                                    format = 'ascii.no_header', delimiter = '\s')

        align_mag_main = Table.read(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[0][0:3]}_named.mag', 
                                    format = 'ascii.no_header', delimiter = '\s')

        align_miss_main = Table.read(f'align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{self.stf_lp_lis[0][0:3]}_named.miss0', 
                                     format = 'ascii.no_header', delimiter = '\s')

        
        # Only want to remove stars in Lp and Kp, and we don't want to subtract calibrator stars
        remove_mask1 = np.array(align_pos_main['col4']) < -10000
        remove_mask2 = np.array(align_pos_main['col2']) < -10000
        remove_mask3 = np.isin(np.array(align_pos_main['col1']).astype(str), self.calibrators)
        remove_stars = (np.array(align_pos_main['col1'])[remove_mask1 + remove_mask2 + remove_mask3]) 
        remove_indices = np.where(np.isin(np.array(align_pos_main['col1']), remove_stars))[0]

        align_pos_main.remove_rows(remove_indices)
        align_mag_main.remove_rows(remove_indices)

        stars = np.array(align_pos_main['col1'])
        x = np.array(align_pos_main['col2'])
        y = np.array(align_pos_main['col3'])
        mag = np.array(align_mag_main['col5'])
        year = np.full(len(mag), align_miss_main['col3'][0])
        flux = np.full(len(mag), 1000000.0)
        fill1 = np.full(len(mag), 28.0)
        fill2 = np.full(len(mag), 1.00)
        fill3 = np.full(len(mag), 5)
        fill4 = np.full(len(mag), 7)
    
        
        input_main = pd.DataFrame({"name": stars, "mag": mag, "year": year, "x": x, "y": y, "fill1": fill1, "fill2": fill2,
                                  "fill3": fill3, "flux": flux, "fill4": fill4})

        input_main.to_csv(f"{os.getcwd()}/force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp_inp.lis", index = False, header = False, sep = '\t')

        
        
        # Run in force mode directory
        os.chdir(f'force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/')
        
        # Main map run
        
        # Image for force mode to work on
        imageFile = f'{os.getcwd()}/mag{self.epoch_lp}_lp.fits'

        stf_inp_lis = f"{os.getcwd()}/mag{self.epoch_lp}_lp_inp.lis"

        fileIDLbatch = f'idlbatch_{self.epoch_lp}_' + 'lp_main'
        fileIDLlog = fileIDLbatch + '.log'
        _batch = open(fileIDLbatch, 'w', )
        _batch.write("find_stf_lis, ")
        _batch.write(f"\"{imageFile}\"" + ", ")
        _batch.write(f"{self.corr}"  + ", ")
        _batch.write(f"\"{stf_inp_lis}\"" + ", ")
        _batch.write(f"deblend={self.deblend}" + ", ")
        _batch.write(f"makePsf={self.makePsf}" + ", ")
        _batch.write(f"makeRes={self.makeRes}" + ", ")
        _batch.write(f"makeStars={self.makeStars}" + ", ")
        _batch.write(f"legacy={self.legacy}" + ", ")
        _batch.write(f"cooStar=\"{self.cooStar}\"" + ", ")
        _batch.write(f"starlist=\"{self.starlist}\"" + ", ")
        _batch.write(f"trimfake={self.trimfake}" + ", ")
        _batch.write(f"force={self.force}" + ", ")
        _batch.write(f"fix_pos={self.fix_pos}")
        _batch.write("\n")
        _batch.write("exit\n")
        _batch.close()
    
        # run forcemode from batch file
        cmd = 'idl < ' + fileIDLbatch + ' >& ' + fileIDLlog   
        print(cmd)
        subprocess.call(cmd, shell=True)




        # Define submap IDs
        submap_ids = [1, 2, 3]
        
        # Loop through each submap
        for submap_id in submap_ids:
            # File suffix
            suffix = self.stf_lp_lis[submap_id][0:5]
        
            # Read tables
            align_pos = Table.read(
                f'{self.working_dir}/align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{suffix}_named.pos',
                format='ascii.no_header', delimiter='\s'
            )
            align_mag = Table.read(
                f'{self.working_dir}/align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{suffix}_named.mag',
                format='ascii.no_header', delimiter='\s'
            )
            align_miss = Table.read(
                f'{self.working_dir}/align/{self.epoch_lp}_v{self.stf_lp_ver}/align_kp_{suffix}_named.miss0',
                format='ascii.no_header', delimiter='\s'
            )
        
            # Identify and remove invalid rows
            remove_mask1 = np.array(align_pos['col4']) < -10000
            remove_mask2 = np.array(align_pos['col2']) < -10000
            remove_mask3 = np.isin(np.array(align_pos['col1']).astype(str), self.calibrators)
            remove_stars = np.array(align_pos['col1'])[remove_mask1 | remove_mask2 | remove_mask3]  # Combine masks with OR
            remove_indices = np.where(np.isin(np.array(align_pos['col1']), remove_stars))[0]
            align_pos.remove_rows(remove_indices)
            align_mag.remove_rows(remove_indices)
        
            # Extract data and create filled arrays
            stars = np.array(align_pos['col1'])
            x = np.array(align_pos['col2'])
            y = np.array(align_pos['col3'])
            mag = np.array(align_mag['col5'])
            year = np.full(len(mag), align_miss['col3'][0])
            flux = np.full(len(mag), 1000000.0)
            fill1 = np.full(len(mag), 28.0)
            fill2 = np.full(len(mag), 1.00)
            fill3 = np.full(len(mag), 5)
            fill4 = np.full(len(mag), 7)
        
            # Create DataFrame
            input_df = pd.DataFrame({
                "name": stars, "mag": mag, "year": year, "x": x, "y": y, 
                "fill1": fill1, "fill2": fill2, "fill3": fill3, "flux": flux, "fill4": fill4
            })
        
            # Save to file
            input_df.to_csv(
                f"{self.working_dir}/force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{submap_id}_inp.lis",
                index=False, header=False, sep='\t'
            )

        

        for i in range(1,3+1):
           # Image for force mode to work on
        
            imageFile = f"{os.getcwd()}/m{self.epoch_lp}_lp_{i}.fits"
            
            stf_inp_lis = f"{os.getcwd()}/m{self.epoch_lp}_lp_{i}_inp.lis"

            fileIDLbatch = f'idlbatch_{self.epoch_lp}_' + f'lp_{i}'
            fileIDLlog = fileIDLbatch + '.log'
            _batch = open(fileIDLbatch, 'w', )
            _batch.write("find_stf_lis, ")
            _batch.write(f"\"{imageFile}\"" + ", ")
            _batch.write(f"{self.corr}"  + ", ")
            _batch.write(f"\"{stf_inp_lis}\"" + ", ")
            _batch.write(f"deblend={self.deblend}" + ", ")
            _batch.write(f"makePsf={self.makePsf}" + ", ")
            _batch.write(f"makeRes={self.makeRes}" + ", ")
            _batch.write(f"makeStars={self.makeStars}" + ", ")
            _batch.write(f"legacy={self.legacy}" + ", ")
            _batch.write(f"cooStar=\"{self.cooStar}\"" + ", ")
            _batch.write(f"starlist=\"{self.starlist}\"" + ", ")
            _batch.write(f"trimfake={self.trimfake}" + ", ")
            _batch.write(f"force={self.force}" + ", ")
            _batch.write(f"fix_pos={self.fix_pos}")
            _batch.write("\n")
            _batch.write("exit\n")
            _batch.close()
            
            # run forcemode from batch file
            cmd = 'idl < ' + fileIDLbatch + ' >& ' + fileIDLlog   
            print(cmd)
            subprocess.call(cmd, shell=True)

        # Save final list of starplanting sources
    def assemble_final_list(self):

        print("""
        
        Assembling final list by only considering area of interest...
        
        """)
        initial_main = Table.read(f"{self.working_dir}/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[0]}_model.lis", delimiter = '\t', format = 'ascii.no_header')
        add_main = Table.read(f"{self.working_dir}/force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp_{self.corr}_stf.lis", delimiter = '\s', format = 'ascii.no_header')
        
        del add_main['col3']
        del add_main['col6']
        del add_main['col7']
        del add_main['col8']

        add_main.rename_column('col4', 'col3')
        add_main.rename_column('col5', 'col4')
        add_main.rename_column('col9', 'col5')

        main_final = vstack([initial_main, add_main])

        mask = (
        (main_final['col3'] >= self.sgra_loc[0] - self.xpixels) & 
        (main_final['col3'] <= self.sgra_loc[0] + self.xpixels) & 
        (main_final['col4'] >= self.sgra_loc[1] - self.ypixels) & 
        (main_final['col4'] <= self.sgra_loc[1] + self.ypixels)
        )
        main_final = main_final[mask]
        

        main_final.write(f"{self.working_dir}/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[0]}_model.lis", delimiter = '\t', format = 'ascii.no_header', overwrite = True)

        for submap_id in range(1, 4):
            # Read the initial and add files for each submap
            initial = Table.read(f"{self.working_dir}/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[submap_id]}_model.lis", delimiter='\t', format='ascii.no_header')
            add = Table.read(f"{self.working_dir}/force_mode/{self.epoch_lp}/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{submap_id}_{self.corr}_stf.lis", delimiter='\s', format='ascii.no_header')
            
            # Remove unnecessary columns from the additional table
            add.remove_column('col3')
            add.remove_column('col6')
            add.remove_column('col7')
            add.remove_column('col8')
    
            # Rename columns for consistency
            add.rename_column('col4', 'col3')
            add.rename_column('col5', 'col4')
            add.rename_column('col9', 'col5')
    
            # Combine the initial and additional submap tables
            combined = vstack([initial, add])

            mask = (
            (combined['col3'] >= self.sgra_loc[0] - self.xpixels) & 
            (combined['col3'] <= self.sgra_loc[0] + self.xpixels) & 
            (combined['col4'] >= self.sgra_loc[1] - self.ypixels) & 
            (combined['col4'] <= self.sgra_loc[1] + self.ypixels)
            )
            combined = combined[mask]
            
            # Save the modified submap table back to the file
            combined.write(f"{self.working_dir}/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[submap_id]}_model.lis", format='ascii.no_header', delimiter='\t', overwrite=True)
                    
    


    def starplant(self):

        print("""
        
        Creating model sources with StarPlanting...
        
        """)
        # Main map
        os.makedirs(f'{self.working_dir}/lp_model/{self.epoch_lp}/', exist_ok = True)
        # image to plant on
        shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/mag{self.epoch_lp}_lp_back.fits', self.working_dir + f'/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[0]}_model.fits')

        # lp psf model
        shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/mag{self.epoch_lp}_lp_psf.fits', self.working_dir + f'/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[0]}_psf.fits')
        
        # Get image size
        with fits.open(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/mag{self.epoch_lp}_lp_back.fits') as hdu:
            hdr = hdu[0].header
            imgsize = hdr['NAXIS1']

        # starplant, d, int_num, rootdir, psfFile, epoch, dr_dir, imsize, filt, stf_lp_lis, input
    
        fileIDLbatch = f'idlbatch_{self.epoch_lp}_main_' + 'lp'
        fileIDLlog = fileIDLbatch + '.log'
        _batch = open(fileIDLbatch, 'w', )
        _batch.write("starplant, ")
        _batch.write('0'  + ", ") # d
        _batch.write('0' + ", ") # int_num
        _batch.write("'" + self.working_dir + f'/lp_model/' + "'" + ', ') # rootdir
        _batch.write("'" + self.working_dir + f'/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[0]}_psf.fits' +"'" + ", ") # psfFile
        _batch.write("'" + f'{self.epoch_lp}' + "'" + ", ") # epoch
        _batch.write("'" + self.working_dir + f'/lp_model/' + "'" + ', ') # dr_dir -- actually has no purpose, remove later
        _batch.write(str(imgsize) + ", ") # imsize
        _batch.write("'" +'lp' + "'" + ", ") # filt
        _batch.write("'" + f'{self.stf_lp_lis[0]}' + "'" + ", ") # stf_lp_lis
        _batch.write("'" + f'{self.working_dir}/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[0]}_model.lis' + "'") # input
        _batch.write("\n")
        _batch.write("exit\n")
        _batch.close()
        
    
        # run starplanting from batch file
        cmd = 'idl < ' + fileIDLbatch + ' >& ' + fileIDLlog   
        subprocess.call(cmd, shell=True)
        if os.path.exists(f'{self.working_dir}/starplanting_batch_files/{fileIDLbatch}'):
            os.remove(f'{self.working_dir}/starplanting_batch_files/{fileIDLbatch}')
            os.remove(f'{self.working_dir}/starplanting_batch_files/{fileIDLlog}')

        shutil.move(f'{self.working_dir}/{fileIDLbatch}', f'{self.working_dir}/starplanting_batch_files/')
        shutil.move(f'{self.working_dir}/{fileIDLlog}', f'{self.working_dir}/starplanting_batch_files/')


    # Sub maps
        for submap_id in range(1, 4):
            # Copy the back and PSF model files for each submap
            shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/m{self.epoch_lp}_lp_{submap_id}_back.fits', 
                        self.working_dir + f'/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[submap_id]}_model.fits')
    
            shutil.copy(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/m{self.epoch_lp}_lp_{submap_id}_psf.fits', 
                        self.working_dir + f'/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[submap_id]}_psf.fits')
    
            # Get image size for the submap
            with fits.open(self.lp_dir + f'starlists/combo/{self.epoch_lp}/starfinder_v{self.stf_lp_ver}/m{self.epoch_lp}_lp_{submap_id}_back.fits') as hdu:
                hdr = hdu[0].header
                imgsize = hdr['NAXIS1']
    
            # Create IDL batch file for each submap
            fileIDLbatch = f'idlbatch_{self.epoch_lp}_{submap_id}_lp'
            fileIDLlog = fileIDLbatch + '.log'
            
            # Write the batch file for running starplanting
            with open(fileIDLbatch, 'w') as _batch:
                _batch.write("starplant, ")
                _batch.write('0' + ", ")
                _batch.write('0' + ", ")
                _batch.write("'" + self.working_dir + f'/lp_model/' + "'" + ', ')
                _batch.write("'" + self.working_dir + f'/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[submap_id]}_psf.fits' + "'" + ", ")
                _batch.write("'" + f'{self.epoch_lp}' + "'" + ", ")
                _batch.write("'" + self.working_dir + f'/lp_model/' + "'" + ', ')
                _batch.write(str(imgsize) + ", ")
                _batch.write("'" + 'lp' + "'" + ", ")
                _batch.write("'" + f'{self.stf_lp_lis[submap_id]}' + "'" + ", ")
                _batch.write("'" + f'{self.working_dir}/lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[submap_id]}_model.lis' + "'")
                _batch.write("\n")
                _batch.write("exit\n")
    
            # Run starplanting from batch file
            cmd = f'idl < {fileIDLbatch} >& {fileIDLlog}'
            subprocess.call(cmd, shell=True)

            if os.path.exists(f'{self.working_dir}/starplanting_batch_files/{fileIDLbatch}'):
                os.remove(f'{self.working_dir}/starplanting_batch_files/{fileIDLbatch}')
                os.remove(f'{self.working_dir}/starplanting_batch_files/{fileIDLlog}')
            shutil.move(f'{self.working_dir}/{fileIDLbatch}', f'{self.working_dir}/starplanting_batch_files/')
            shutil.move(f'{self.working_dir}/{fileIDLlog}', f'{self.working_dir}/starplanting_batch_files/')
                
                
        
                
    def create_residual(self):   

        print("""
        Creating residual and saving them to combo directories...
        
        """)
        os.chdir(self.working_dir)
        os.makedirs(self.working_dir + f'lp_residuals/{self.epoch_lp}/', exist_ok = True)
        os.makedirs(self.working_dir + f'lp_combo/{self.epoch_lp}/', exist_ok = True)
        # copy Lp combo
        shutil.copy(self.lp_dir + f'/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp.fits', self.working_dir + f'/lp_combo/{self.epoch_lp}/{self.epoch_lp}_combo.fits')
         
        # Lp combo
        data_combo, hdr = fits.getdata(f'lp_combo/{self.epoch_lp}/{self.epoch_lp}_combo.fits', header = True)
            
        # Lp model
        data_model = fits.getdata(f'lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[0]}_model.fits')
            
        # Lp residual
        data_residual = data_combo - data_model
        
        
        # Write fits file
        fits.writeto(f'lp_residuals/{self.epoch_lp}/{self.epoch_lp}_v{self.stf_lp_ver}_{self.stf_lp_lis[0]}_residual.fits', data_residual, hdr, overwrite = True)

        shutil.copy(f'lp_residuals/{self.epoch_lp}/{self.epoch_lp}_v{self.stf_lp_ver}_{self.stf_lp_lis[0]}_residual.fits', 
                   self.lp_dir + f'/combo/{self.epoch_lp}/mag{self.epoch_lp}_lp.fits')



        for i in range(1,4):
            shutil.copy(self.lp_dir + f'/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}.fits', self.working_dir + f'/lp_combo/{self.epoch_lp}/{self.epoch_lp}_{i}.fits')
             
            # Lp combo
            data_combo, hdr = fits.getdata(f'lp_combo/{self.epoch_lp}/{self.epoch_lp}_{i}.fits', header = True)
                
            # Lp model
            data_model = fits.getdata(f'lp_model/{self.epoch_lp}/{self.epoch_lp}_{self.stf_lp_lis[i]}_model.fits')
                
            # Lp residual
            data_residual = data_combo - data_model
            
            
            # Write fits file
            fits.writeto(f'lp_residuals/{self.epoch_lp}/{self.epoch_lp}_v{self.stf_lp_ver}_{self.stf_lp_lis[i]}_residual.fits', data_residual, hdr, overwrite = True)

            shutil.copy(f'lp_residuals/{self.epoch_lp}/{self.epoch_lp}_v{self.stf_lp_ver}_{self.stf_lp_lis[i]}_residual.fits', 
                        self.lp_dir + f'/combo/{self.epoch_lp}/m{self.epoch_lp}_lp_{i}.fits')

            
                    
                        
            
                        
                        
                        
