# G2 Evolution Project
The purpose of this project is to investigate the evolution of various properties of G2. 
The Astrometry is measured with NIRC2 (Lp ~ 3.8µm) and OSIRIS (Br_gamma_ ~ 2.166µm) from the years 2002 - 2023 (NIRC2)
and 2006 - 2024 (OSIRIS). We also measure its radial velocity with OSIRIS and photometry with NIRC2. 




## Lp residual analysis code
To obtain NIRC2 astrometric and photometric measurements, we employ the **analysis_residuals.py** object.
The **Lp_residual_analysis.py** code utilizes the analysis_residuals.py object.

Purpose: The purpose of this code is to un-confuse G2 with Kp sources in Lp images and then to subsequently 
obtain its astrometry and photometry with starfinder. 


The outline of this code is as follows:

Create directory to work in, this will be referenced as the working directory (e.g. `/g/ghez/seanz/projects/G-objects/starplanter/`)

1. Setup align directory (`def setup_align(self)`):
	1. make align directory `align/{epoch_lp}_v{version}/`
	2. copy `orbits.dat` and `label.dat` into align directory
	3. Write `align_kp_{corr}_named.lis` file with main map calibration list written inside for `lp` and `kp`
	4. copy main map cal.lis files (`/g/ghez/data/dr/dr2/starlists/combo/{epoch}/starfinder_v{version}/mag{epoch}_kp_0.8_cal.lis`)into align directory
	5. Repeat steps 2-4 for submaps
2. Run align between `lp` and `kp` starlists (`def run_align(self)`):
	1. change directory to the align directory 
	2. Run java align for main map and submaps
3. Setup force mode directory (`def setup_force(self)`):
	1. change directory back to working directory
	2. make force mode directory `force_mode/{epoch_lp}/combo/{epoch_lp}/`
	3. copy over main map (`mag{epoch_lp}_lp.fits`) into force directory
	4. copy over the input list for force mode which will be the missed kp stars in the original align. This is done so that we are more robust in the removal of kp sources in lp image (`align/{epoch_lp}_v{version}/align_kp_0.7_named.miss0`)
	5. copy over the necessary files for force mode to run (`back.fits`, `psf.fits`, `.coo`)
	6. repeat steps 3-5 for submaps
4. run force mode (`def missed_kp(self)):
	1. Run the idl procedure `find_stf_lis` on the main and sub maps in lp. 
	2. Have force = true, and the input list be specified as the align.miss file we copied over in step 4 in the setup force directory section. Rename it to `mag{epoch_lp}_lp_inp.lis` before running.
5. Generate initial starplanting list by inputing the missed kp sources that are now detected in lp image (`def align_force_compare(self)`)
	1. At this stage, two sets of starlists have been generated: `align.pos` files and `forced_stf.lis` files. We want to create a robust list of stars with positions, calibrated magnitudes and fluxes to starplant in lp. 
	2. `def align_force_compare(self)` will find stars that are matched between `align.pos` and `forced_stf.lis` files for main and sub maps. We want to have only the sources detected in `forced_stf.lis` to be our initial starplanting list
	3. find matches between `forced_stf.lis` and `align.pos`, remove those stars from `forced_stf.lis`
	4. save edited `forced_stf.lis` as initial model list --> `{epoch_lp}_0.7_stf_model.lis` with name, magnitude, x, y, and flux as columns.
	5. Repeat steps 3&4 for submaps
6. Get magnitudes and fluxes for kp/lp sources from the lp image `def kp_and_lp_sources(self):`
	1. from the `align.pos` file, remove stars that are not matched in both lp and kp
	2. with these stars, find the magnitudes from the `align.mag` file
	3. create starfinder force mode input list that contains the star names, x, y, magnitude, year, flux (any float, starfinder will find this), and 4 random filled columns containing float integers
	4. Save this as a starfinder force mode input list
	5. Run starfinder force mode
	6. repeat steps 1-5 for submaps
7. Construct final starplanting list `def assemble_final_list(self):`
	1. combine the initial starplanting list and the new starfinder force mode output list for main and submaps
8. Starplant a model image `def starplant(self):`
	1. make lp_model directory --> `{working directory}/lp_model/{epoch_lp}`
	2. copy over the corresponding background fits image 
	3. copy over the corresponding psf model fits image
	4. obtain the corresponding image size
	5. run starplanting 
	6. repeat steps 2-5 for submaps
9. Create residuals 'def create_residual(self):'
	1. Copy over Lp combo image
	2. Subtract combo and starplanted iamge
	3. Replace resulting image in data combo directory
	4. Repeat steps 1-3 for submaps.
	
	
	





## Maintainers
- Sean Granados: 
	Affiliation: UCLA
	Email: (seangranados@astro.ucla.edu)


## Last updated
- 11/27/2024 -- Sean Granados: Created document
