FUNCTION psf_clean, image, psf, x, y, fitbox, _EXTRA = extra
END

function find_stf_lis_input_list, input
    ; Count the number of stars to deal with
    numStars= file_lines(input)

    ; If something unexpected happens get out
    if (numStars lt 1) then begin
        print, "****** 0 stars in the input list. SOMETHING WRONG. ******"
        retall
     endif
    
    ; Create variables to store stuff
    inName = strarr(numStars)
    inMag = fltarr(numStars)
    inDate = fltarr(numStars)
    inX = fltarr(numStars)
    inY = fltarr(numStars)
    inCorr = fltarr(numStars)

    ;Read in the File
    openr, u, input, /get_lun

    ; Loop through the file and read in all the info
    for i=0, numStars-1 do begin
        temp = ''
        readf, u, temp
        line = strsplit(temp, /EXTRACT)
        numstrings = N_ELEMENTS(line)
        ; Two possible formats in reading, with or without error
        ; and with or without the number at the end

        ; These are the same in either case
        inName[i]= string(line[0])
        inMag[i] = float(line[1])
        inDate[i] = float(line[2])
        inX[i] = float(line[3])
        inY[i] = float(line[4])

        ; Depending on if the file includes errors or not
        ; pick out the correct correlation values
        if (numstrings eq 11 or numstrings eq 12) then begin
            inCorr[i] = float(line[8])
         endif else if (numstrings eq 10 or numstrings eq 9 ) then begin
            inCorr[i] = float(line[6])
         endif else begin
            ; If the input file is an unrecognizable format, stop
            ; if (numstrings lt 10 or numstrings gt 12) then begin
            print, "Input file is an unrecognized format - retalling."
            retall
         endelse

      endfor
    close, u & free_lun, u

    ; Sort by Magnitude BUT keep the first star at the top of the list
    order = sort(inMag[1:*])
    order = [0, order]

    inStars = REPLICATE( { name: '', $
            mag: 0.0, $
            x: 0.0, $
            y: 0.0, $
            c: 0.0, $
            year: 0.0, $
            i0: 0.0 $
          }, N_ELEMENTS(inX))
    instars.name = inName[order]
    instars.x = inX[order]
    instars.y = inY[order]
    instars.c = inCorr[order]
    instars.year = inDate[order]
    instars.mag = inMag[order]

    return, instars
end

;+
; NAME:
;    sigFile = get_sig_file(inRoot)
;
; PURPOSE:
;    Get the name of the significance file or "single"
;-
FUNCTION get_sig_file, inRoot
    sigFile = inRoot + '_sig.fits'
    if (exist(sigFile) EQ 0) then sigFile = inRoot + '_wgt.fits'
    if (exist(sigFile) EQ 0) then sigFile = inRoot + '.sig'
    if (exist(sigFile) EQ 0) then sigFile = 'single'

    return, sigFile
END

;+
; NAME:
;    maxThreshold = get_max_threshold(maxFile)
;
; PURPOSE:
;    Get the saturation limit from a *.max file if it exists.
;-
FUNCTION get_max_threshold, image, maxFile
    if exist(maxFile) then begin
        openr, _max, maxFile, /get_lun
        readf, _max, maxThreshold
        free_lun, _max
        print, format='(%"FIND_STF: Saturation above %d")', maxThreshold
    endif else begin
        print, format='(%"FIND_STF: Assuming no saturation issues")'
        ;; There was a super-sneaky bug here... the max threshold must
        ;; NEVER be crossed even if weird, negative backgrounds are
        ;; found by the fitter. Otherwise, you get a hole in the PSF.
        ;; I am setting it way high just in case. JLU 2018-10-30
        maxThreshold = max(image) * 2.0
    endelse

    return, maxThreshold
END

;+
; NAME:
;    psf_size = get_psf_size(psf_size_keyword, default_size, scale)
;
; PURPOSE:
;    Get the psf size if it is set. Otherwise, use the default.
;    PSF size is in arcsec.
;-
FUNCTION get_psf_size, psf_size_keyword, default_size, scale
    ; Setup the PSF size
    if (psf_size_keyword EQ 0) then begin
       psf_size = default_size           ; arcsec
    endif else begin
       psf_size = psf_size_keyword
    endelse
    psf_size = floor(psf_size / scale)
    print, 'FIND_STF: psf_size = ', psf_size

    return, psf_size
END

;+
; NAME:
;    noFixPsf = get_pix_sf(isSpeckle, scale)
;
; PURPOSE:
;    Determine whether we will do PSF halo smoothing and PSF fixing.
;
; RETURN:
;    noFixPsf: 1 if we will not be fixing the PSF. 
;-
FUNCTION get_fix_psf, isSpeckle, scale
    if (isSpeckle EQ 1 or scale GT 0.015) then begin
        no_psf_fix = 0            ; speckle or seeing data needs halo smoothing
    endif else begin
        no_psf_fix = 1
    endelse

    return, no_psf_fix
END


;+
; NAME:
;    find_stf_lis, image, corr, input, year=year, 
;              /deblend, /isSpeckle,
;              /makePsf, /trimframes, 
;              [cooStar=starName]
;
; PURPOSE:
;    Run starfinder on the specified image (FITS format) using a
;    previously created PSF (FITS format) and background (FITS). 
;    Specify the correlation
;    value and whether to turn deblending on.
;
; INPUTS:
;    image -- The name of the fits file for the image to work on. For
;             an image called 03jun/combo/mag03jun.fits, the
;             PSF and background images should be named
;             03jun/combo/mag03jun_psf.fits and
;             03jun/combo/mag03jun_back.fits.  
;    corr -- An array of correlation values. Starfinder will be run
;            once at each of the correlation values specified.
;    input -- Input list of stars.   
; OPTIONAL INPUTS:
;    year -- Specify the date in fractional years (with at least 3 
;            decimal places). If not specified, then the date will be
;            determined from the header of the input image.
;    /deblend -- Default value is to turn deblending off. 
;    /isSpeckle -- Will use *.sig file if turned on or will look for 
;                  a *_sig.fits file otherwise. If there is no *_sig.fits
;                  file then check for a *_wgt.fits (drizzle output).
;                  Default is 0.
;    /makePsf -- Will make a new PSF from several named sources within
;                the field (see stf_psf_positions.pro for
;                details). WARNING: This will write over the previously
;                existing PSF. This will also make a new background
;                for the image. WARNING: This will write over the 
;                previously existing backgrounds.
;                If aoopt = 1 and makePsf = 0 then makeGrid = 0 (see below).
;                Default is makePsf = 0.
;    /makeRes -- Make a brand new residual image (image - stars - background). The
;                file is stored in the same directory as the image to be
;                worked on and has a file name such as mag04jul_res.fits
;                Default is 0.
;    /makeStars -- Make a brand new stars image. The
;                file is stored in the same directory as the image to be
;                worked on and has a file name such as
;                mag04jul_stars.fits
;                Default is 0.
;    /quick -- Use this for a quick and dirty star detection suitable
;              for aligning single frames in order to get the proper
;              offsets for combining. Astrometry and photometery of
;              individual stars and PSFs should not be trusted. This
;              sets a very high sigma threshold for star detection and
;              does the minimal number of iterations even to create a
;              PSF. Only ever use quick on non-speckle data.
;     /trimframes -- Use this to trim the edges of stars without
;                    enough frames in them. Uses a flat drop flat 
;                    algorithm to keep high correlation sources
;                    with less frames than low correlation sources 
;                    at the same frames
;    cooStar -- The name of the star whose coordinates are located
;               in the *.coo file. If not specified, then the star
;               is assumed to be "irs16C". This is neccesary in order
;               to do automated PSF extraction with known stars. See
;               stf_psf_positions to learn more about this.
;               Default is 'irs16C'.
;    /trimfake -- Use this to measure the PSF envelope from plots of 
;                 pairwise delta-magnitudes vs. separation. Once the 
;                 envelope is determined, remove all sources that fall
;                 too close and too faint to existing sources. 
;                 (Default = 0, no trim)
;    backBoxFWHM -- Box size (in units of PSF FWHM) for the
;                   background.
;    psfSearchBox -- Box size (in arcsec) for centroiding the PSFs in
;              stf_psf_positions.pro
;    secSearchBox -- Box size (in arcsec) for centroiding the PSFs in 
;              stf_psf_positions.pro
;    psfSize -- Set psf size in arcsec; if 0 the default PSF size is used 
;               (3 arcsec, 301 pixels). 
;               Default is 0.
;    /fixPsf -- Set to 1 for standard PSF postprocessing (circular
;               mask, radial and azimuthal smoothing).
;               Default is 0.
;    starlist -- PSF reference source list used by stf_psf_positions.
;                Default is '/g/ghez/data/gc/source_list/psf_central.dat'.
;    backBoxFWHM -- Spatial scale of the background estimation in 
;                   psf_extract and starfinder in units of FWHM of PSF.
;                   Default is 9.
;    /legacy -- Use this to run AIROPA in legacy mode (Starfinder
;               1.6). This mode observes all of the keywords above,
;               but uses the earlier default settings for all keywords
;               (see header of find_stf_legacy). If legacy = 0 and
;               aoopt = 0 then AIROPA is executed in 'single' PSF
;               AIROPA mode (different PSF extraction and postprocessing).
;    /aoopt -- Use this to run ARIOPA in variable PSF mode; many of
;              the followingkeywords set different configurations for this mode;
;              Default is aoopt = 0
;    part_size -- applies only if aoopt = 1; step width of PSF
;                 grid; used in sf_partition.pro to generate PSF grid
;                 and grid positions according to Starfinder
;                 conventions. part_size must be set if aoopt = 1 and
;                 partition = 0.
;    /makeGrid -- applies only if aoopt = 1;
;            Set to 1 to generate the OTF ratio grid from instrumental
;            OTF maps and MASS DIMM files.
;            Set to 0 if OTF grid of previous run is available; then 
;            grid generation is skipped and files following standard 
;            naming schemes are read or input from:
;                 otfgridFile_real, otfgridFile_imag, posFile
;            If aoopt = 1 and makePsf = 0 then makeGrid is
;            automatically set to 0.
;   inst_grid_name_real -- applies only if aoopt = 1; provide instrumental
;                     OTF grid name of the real part if generated independently.
;                     Default is 'undefined'
;   inst_grid_name_imag -- applies only if aoopt = 1; provide instrumental
;                     OTF grid name of the imaginary part if generated independently.
;                     Default is 'undefined'
;   corr_file -- applies only if aoopt = 1; only required if makeGrid = 1,
;                and otfgridFile_real and otfgridFile_imag undefined;
;                provides path and name to correction file to correct
;                for low-level differences between atm. OTF generated
;                at 301 pixels vs. atm. OTF generated at 101 pixels
;                with resampling to 301.
;   dimm_file -- applies only if aoopt = 1; only required if makeGrid = 1,
;                and otfgridFile_real and otfgridFile_imag undefined;
;                path and name of DIMM file
;   mass_file -- applies only if aoopt = 1; only required if makeGrid = 1,
;                and otfgridFile_real and otfgridFile_imag undefined;
;                path and name of MASS file
;   ttStar -- The name of the TT star used in aoopt=1 mode. The star
;             must be located in the psf starlist. 
;   gsStar -- The name of the guide star star used in aoopt=1 mode. The star
;             must be located in the psf starlist. If the guide star
;             is not specified and the image was obtained in LGS mode,
;             then, the LGS will be assumed to be at the center of the image.
;   config_file -- applies only if aoopt = 1; only required if makePsf =
;                  1, and otfgridFile_real and otfgridFile_imag
;                  undefined; path and name of the instrumental
;                  configuration file
;   nearestneighbor -- applies only if aoopt = 1; switch to
;                      alternative mechanism to relate star position
;                      to PSF in gridded SF
;                      based on nearest neighbor; default = 0
;   partition -- applies only if aoopt = 1; applies only if
;                nearestneighbor = 1; path and name to partition fits file if
;                generated outside find_stf (faster)
;   n_grid_over -- factor for generating finer grid; interpolation for
;                  atm. OTF ratio is executed pixelwise onto a grid
;                  that has n_grid_over grid points for each original
;                  grid point.
;   varPSF -- applies only if aoopt = 1; experimental!!! disregard for now.
;   otfgridFile_real -- applies only if aoopt = 1 and makeGrid = 0;
;                       otfgridFile_real, _imag, and posFile need to
;                       be set together. Provide externally generated
;                       OTF ratio ; not atm. or instr. OTF ratio
;                       generation is executed.
;   otfgridFile_imag -- applies only if aoopt = 1 and makeGrid = 0;
;                       otfgridFile_real, _imag, and posFile need to be set
;                       together. Provide externally generated OTF ratio; not atm. or
;                       instr. OTF ratio generation is executed.
;   phase_map_folder -- Folder where are stored the phase map files listed in the airopa config file
;   posFile -- applies only if aoopt = 1 and makeGrid = 0; externally generated file with grid positions
;   back_clip -- PSF post-processing. If 1 then noise-dominated areas of the PSF are clipped to zero; default = 0
;   field -- switch to different PSF clipping parameters; experimental
;   weighted -- PSF reference star weighted averaging; default = 0.
;
;   force -- execute StarFinder in forced mode; requires lists of x-y-positions and fluxes; default = 0
;   x_force -- applies only if force = 1; list of x-positions in forced mode (1D array)
;   y_force -- applies only if force = 1; list of y-positions in forced mode (1D array)
;   f_force -- applies only if force = 1; list of fluxes in forced mode (1D array)
;
;   center_by_fitting -- added to _Extra; passed on to
;                        superpose_stars; registeres the PSF reference stars by
;                        fitting each reference star to the first reference star
;   clean_rad -- inherited through _Extra, passed on to superpose_stars;
;                masking pixels with intensity > 0.1 * max(PSF)
;                outside a radius of size clean_rad
;                during PSF re-extraction, default is 'undefined';
;                proposal value is 25.
;   flat -- inherited through _Extra, passed on to
;           estimate_background; force background to be a constant
;           number over the field; default = 0
;   subtract -- inherited through _Extra, passed in to cont_psf;
;               subtracting threshold from clipped PSF before
;               normalization
;   debug -- keyword to run in debug mode; if debug = 1 then
;            intermediated results of the PSF extraction are written
;            out; debug = 2 for experts only
;
; OUTPUT:
;    <img_file>_<corr>_stf.lis
;        A starlist with positions, magnitudes, fluxes, etc.
;
;    <img_file>_<corr>_metrics.txt
;        A text file with the same length as the above starlist. But
;        this file contains the sqrt(FVU) as calculated during the
;        actual PSF fitting process (i.e. including local linear
;        background, etc.). 
; 
; EXAMPLE:
;    find_stf, 'mag04jul.fits', [0.7, 0.6, 0.5], /deblend
;    
;
; MODIFICATION HISTORY:
;    08/05/2004 - Jessica Lu adapted from Seth Hornstein's find_mwsaa.
;    07/07/2005 - Jessica Lu modified to use exsiting PSF and
;                 background images in the same directory as the image
;                 itself.
;    08/10/2005 - Added makeRes and makeStars and saturated flags.
;    09/19/2005 - Added trimframes (M. Rafelski)
;    09/17/2007 - Added a flag to control which PSF starlist is read in to
;                 stf_psf_positions.pro. This is helpful for Arches 
;                 and M92 data. Default is set to read in the GC PSF 
;                 positions (psfstars/psf_central.dat). (S. Yelda)
;    08/20/2009 - See SVN for modifications history.
;
;-
pro find_stf_lis, imageFile, corr, input, $
              deblend=deblend, year=year, isSpeckle=isSpeckle, quick=quick, $
              makePsf=makePsf, makeRes=makeRes, makeStars=makeStars, $
              makeRep=makeRep, legacy = legacy, trimframes=trimframes, trimfake=trimfake, $
              cooStar=cooStar, guideStar=guideStar, ttStar=ttStar, $
              starlist=starlist, psfSize=psfSize, fixPsf=fixPsf, $
              backBoxFWHM=backBoxFWHM, psfSearchBox=psfSearchBox, secSearchBox=secSearchBox, $
              aoopt = aoopt, part_size = part_size, makeGrid = makeGrid, inst_grid_name_real=inst_grid_name_real, $
              inst_grid_name_imag=inst_grid_name_imag, corr_file = corr_file, dimm_file = dimm_file, mass_file = mass_file, $
              config_file = config_file, center_by_fitting = center_by_fitting, $
              back_clip = back_clip, weighted = weighted, $
              field = field, force = force, x_force = x_force, y_force = y_force, f_force = f_force, add_force = add_force, $
              nearestneighbor = nearestneighbor, partition = partition, varPSF = varPSF, $ 
              otfgridFile_real = otfgridFile_real, otfgridFile_imag = otfgridFile_imag, posFile = posFile, $
              n_grid_over = n_grid_over, phase_map_folder=phase_map_folder, save_otf=save_otf, $
              fix_Psf_Cos = fix_Psf_Cos, fix_Psf_Smooth = fix_Psf_Smooth, fix_Psf_HaloClip = fix_Psf_HaloClip, $
              fix_Psf_Trim = fix_Psf_Trim, fix_Psf_maskrad = fix_Psf_maskrad, fix_Psf_nSigmaStart = fix_Psf_nSigmaStart, $
              _EXTRA = extra, fix_pos = fix_pos

    RESOLVE_ROUTINE, 'starfinder', /compile_full_file
    RESOLVE_ROUTINE, 'find_stf_legacy', /compile_full_file 
    RESOLVE_ROUTINE, 'psf_extract', /compile_full_file
    
    if n_elements(imageFile) EQ 0 then begin
        print, 'FIND_STF: Image file not specified!'
        return
    endif
    
    if n_elements(corr) EQ 0 then begin
        print, 'FIND_STF: Correlation not specified!'
        return
    endif

    if n_elements(input) EQ 0 then begin
        print, 'FIND_STF: Input starlist not specified!'
        return
    endif


    if not (keyword_set(deblend)) then deblend=0
    if not (keyword_set(makePsf)) then makePsf=0
    if not (keyword_set(makeGrid)) then makeGrid=0
    if not (keyword_set(makeRes)) then makeRes=0
    if not (keyword_set(makeStars)) then makeStars=0
    if not (keyword_set(makeRep)) then makeRep=0
    if not (keyword_set(trimfake)) then trimfake=0
    if not (keyword_set(aoopt)) then aoopt=0
    if not (keyword_set(legacy)) then legacy=0
    if not (keyword_set(isSpeckle)) then isSpeckle=0
    if not (keyword_set(cooStar)) then cooStar='irs16C'
    if not (keyword_set(psfSize)) then psfSize = 0
    if not (keyword_set(fixPsf)) then fixPsf = 0
    if not keyword_set(starlist) then starlist = '/g/ghez/data/gc/source_list/psf_central.dat'
    if not keyword_set(backBoxFWHM) then backBoxFWHM = 9
    if (aoopt eq 1 and makePsf eq 0) then makeGrid = 0
    if not (keyword_set(back_clip)) then back_clip = 0
    if not (keyword_set(weighted)) then weighted = 0
    if not (keyword_set(field)) then field = 0
    if not (n_elements(save_otf)) then save_otf = 1
    if not (keyword_set(n_grid_over)) then n_grid_over = 5L
    if not (keyword_set(force)) then force = 0
    if not (keyword_set(fix_pos)) then fix_pos = 0
;    if not (keyword_set(force)) then begin
;       force = 0
;    endif else begin
;       if (n_elements(x_force) eq 0 OR n_elements(y_force) eq 0 OR n_elements(f_force) eq 0) then begin
;          print, 'FIND_STF: For forced modes please provide [x,y] position and fluxes'
;          return
;       endif else begin
;          x = x_force
;          y = y_force
;          f = f_force
;       endelse
;    endelse    
    if not (keyword_set(add_force)) then begin
       add_force = 0
    endif else begin
       if (n_elements(x_force) eq 0 OR n_elements(y_force) eq 0 OR n_elements(f_force) eq 0) then begin
          print, 'FIND_STF: For forced modes please provide [x,y] position and fluxes'
          return
       endif 
    endelse
    
    print, 'FIND_STF: Starting at ', SYSTIME()
    print, 'FIND_STF: Using PSF starlist ', starlist
    print, 'FIND_STF: trimfake = ', trimfake
    print, 'FIND_STF: aoopt = ', aoopt
;STOP
    if legacy eq 1 then begin
       find_stf_legacy, imageFile, corr, $
                        deblend=deblend, year=year, isSpeckle=isSpeckle, quick=quick, $
                        makePsf=makePsf, makeRes=makeRes, makeStars=makeStars, makeRep=makeRep, $
                        trimframes=trimframes, trimfake=trimfake, $
                        cooStar=cooStar, starlist=starlist, psfSize=psfSize, fixPsf=fixPsf, $
                        backBoxFWHM=backBoxFWHM, psfSearchBox=psfSearchBox, secSearchBox=secSearchBox, $
                        _EXTRA=extra, fix_pos = fix_pos, force=force, input=input
       return
       
    endif

    ;;------------------------------
    ;;
    ;; START new AIROPA versions
    ;;
    ;;------------------------------

    ;; Print out the files we are working on.
    print, format='(%"FIND_STF: Image      %s")', imageFile
    
    ;; Construct root file name
    startIdx = strpos(imageFile, '/', /reverse_search) + 1
    stopIdx = strpos(imageFile, '.fits', /reverse_search)

    outRoot = strmid(imageFile, startIdx, stopIdx - startIdx)
    inRoot = strmid(imageFile, 0, stopIdx)
    inDir = file_dirname(imageFile, /mark_directory)

    ;; Here are some of the file names.
    bkgFile = inRoot + '_back.fits'
    resFile = inRoot + '_res.fits'
    starsFile = inRoot + '_stars.fits'
    repFile = inRoot + '_rep.fits'
    maxFile = inRoot + '.max'
    cooFile = inRoot + '.coo'

    ;;----------
    ;; Get image and image information
    ;;----------
    fits_read, imageFile, image, hdr, /pdu
    siz = size52(image, /dim)

    ;; Get image parameters
    theta = get_pa(hdr)            ; good for HST, NIRC2, speckle, OSIRIS
    scale = get_scale(hdr)         ; arcsec per pixel
    filter = get_filter_name(hdr)  ; should match in psf_central.dat
    year = calc_year(hdr)          ; get the decimal year the image was taken
    date_obs = repstr( strcompress( sxpar(hdr, 'DATE-OBS'), /remove_all), "-", "") ; e.g. 20180723
    psf_size = get_psf_size(psfSize, 2.99, scale)    ; Setup the PSF size
    maxThreshold = get_max_threshold(image, maxFile) ; Get max DN before saturation
    refCoo = get_coo_star_xy(cooFile, cooStar)       ; pixel coorindates of the COO star
    
    print, format='(%"FIND_STF: PA = %6.1f, FILTER = %s, SCALE = %5.3f, YEAR = %8.3f")', $
           theta, filter, scale, year

    if (keyword_set(aoopt)) then begin
       if n_elements(partition) EQ 0 AND n_elements(part_size) EQ 0 then begin
          print, 'FIND_STF: Partition size not specified!'
          return
       endif
   
       ;; Get pixel position of guide star (GS) and tip-tilt (TT) star. 
       get_guide_star_xy, starlist, year, scale, theta, siz, refCoo, $
                          cooStar, gsStar, ttStar, $ 
                          gs_x, gs_y, tt_x, tt_y
    
       psfFile = inRoot + '_on_axis_psf.fits'
       psfgridFile = inRoot + '_psf_grid.fits'
       if not keyword_set(otfgridFile_real) then otfgridFile_real = inRoot + '_otf_ratio_grid_real.fits'
       if not keyword_set(otfgridFile_imag) then otfgridFile_imag  = inRoot + '_otf_ratio_grid_imag.fits'
       if not keyword_set(posFile) then gridposFile = inRoot + '_grid_pos.fits' else gridposFile = posFile
       if not keyword_set(mass_file) then mass_file = inDir + date_obs + '.masspro.dat'
       if not keyword_set(dimm_file) then dimm_file = inDir + date_obs + '.dimm.dat'
       if not keyword_set(config_file) then config_file = getenv('AIROPA_DATA_PATH') + '/ref_files/airopa.config'

       ;; In gridded mode (and when we aren't inputting a grid),
       ;; load up the configuration files, OTF ratio grid, and MASS-DIMM data. 
       if keyword_set(makeGrid) then begin
          COMMON share_config, config_structure
          if size(config_structure, /type) ne 8 then begin
              print, 'FIND_STF: Loading config_file = ', config_file
              airopa_config, config_file, path=phase_map_folder
          endif
       endif
    endif else begin
       psfFile = inRoot + '_psf.fits'
    endelse

    ;; Print out the files we are working on.
    if (makePsf EQ 0) then begin
        print, format='(%"FIND_STF: PSF        %s")', psfFile
        print, format='(%"FIND_STF: Background %s")', bkgFile
    endif

    ;; Get the significance file, either size of image, or "single"
    sigFile = get_sig_file(inRoot)

    ;; Report some other values 
    print, format='(%"FIND_STF: Fix PSF = %d")', fixPsf

    ;; Determine if we should smooth the PSF (cos filter + radial
    ;; smoothing) or background clip the PSF in the second round of
    ;; PSF extraction and secondary removal. 
    ;; These are hard-coded right now.
    if (isSpeckle EQ 1) then begin
       ; Speckle case
       noPsfCos = 0
       noPsfSmooth = 0
       noPsfHaloClip = 0
       noPsfTrim = 1
    endif else if (scale GT 0.03) then begin
       ; Seeing-limited case
       noPsfCos = 1
       noPsfSmooth = 0
       noPsfHaloClip = 1
       noPsfTrim = 1
    endif else begin
       ; AO case
       noPsfCos = 0
       noPsfSmooth = 1
       noPsfHaloClip = 0
       noPsfTrim = 1
    endelse

    ;; Set parameters for PSF noise cut and contiguous signal extraction
    ;;
    ;; maskrad: as large as the PSF halo could presumably reach
    ;; n_sigma_1: use only pixels of n_sigma significance.
    ;;    This only has an effect if the maskrad keyword is set.
    ;; If you use subsequently halo smoothing
    ;; then n_sigma can be given a rather low level (0.1/0.5)
    maskrad_1 = floor(0.65 / scale)  ; Adopt 0.65'' as the masking radius.
    n_sigma_1 = 1.0 ; Gunther's value 1.0

    ;; WHY DO WE HAVE A SECOND SET???
    if field gt 0 then begin
        maskrad = floor(1.0 / scale)
    endif else begin
        maskrad = floor(1.05 / scale)
    endelse
    n_sigma_start = 1.0         ; Gunther's value 1.0

    ;; Override fixPSF Parameters at user request.
    ;; Note these only get applied in the 2nd - last iteration (not
    ;;                                          the first).
    if (keyword_set(fix_Psf_Cos)) then noPsfCos = 0 
    if (keyword_set(fix_Psf_Smooth)) then noPsfSmooth = 0 
    if (keyword_set(fix_Psf_HaloClip)) then noPsfHaloClip = 0 
    if (keyword_set(fix_Psf_Trim)) then noPsfTrim = 0
    if (keyword_set(fix_Psf_maskrad)) then maskrad = fix_Psf_maskrad
    if (keyword_set(fix_Psf_nSigmaStart)) then n_sigma_start = fix_Psf_nSigmaStart
    
    
    ;; Set PSF extraction to use median not average. Also increase the
    ;; interpolation and centroiding tolerances.
    if (N_ELEMENTS(extra) NE 0) then begin
        extra = CREATE_STRUCT('AVGTYPE', 2, $
                              'POS_TOL', 0.0001, $
                              'ASTROMETRIC_TOL', 0.0001, $
                              'INTERP_TYPE', 'SPLINE3', $
                              extra)
    endif else begin
        extra = CREATE_STRUCT('AVGTYPE', 2, $
                              'POS_TOL', 0.0001, $
                              'ASTROMETRIC_TOL', 0.0001, $
                              'INTERP_TYPE', 'SPLINE3')
    endelse
    
    ; Set parameter for PSF interpolation
    grid_over = n_grid_over

    if (keyword_set(aoopt)) then begin
       UTC_date = get_utc_date(hdr)
       cenwave = get_lambda(hdr)
       parang = get_parang(hdr)
       lgs_alt = get_lgs_alt(hdr)
       zen = get_zen(hdr)
       
       if keyword_set(makeGrid) then begin
          if config_structure.ao_mode eq 'NGS' then begin
             grid_input_structure = {   image_size:    siz[0], $
                                        psf_size:      psf_size, $
                                        grid_size:     part_size, $
                                        grid_over:     grid_over, $
                                        lambda:        cenwave, $
                                        rotposn:       theta, $
                                        parang:        parang, $
                                        zen:           zen, $
                                        lgsalt:        lgs_alt, $
                                        date:          UTC_date, $
                                        gs_x:          gs_x, $
                                        gs_y:          gs_y, $
                                        mass_file:     mass_file, $
                                        dimm_file:     dimm_file}
         
          endif else begin
             grid_input_structure = {   image_size:    siz[0], $
                                        psf_size:      psf_size, $
                                        grid_size:     part_size, $
                                        grid_over:     grid_over, $
                                        lambda:        cenwave, $
                                        rotposn:       theta, $
                                        parang:        parang, $
                                        zen:           zen, $
                                        lgsalt:        lgs_alt, $
                                        date:          UTC_date, $
                                        ttgs_x:        tt_x, $
                                        ttgs_y:        tt_y, $
                                        lgs_x:         gs_x, $
                                        lgs_y:         gs_y, $
                                        mass_file:     mass_file, $
                                        dimm_file:     dimm_file}
             atm_corr = readfits(corr_file)
          endelse
   
          if keyword_set(inst_grid_name_real) then begin
              print, 'FIND_STF: Reading instrumental OTF grid as input. '
              print, '   ', inst_grid_name_real
              print, '   ', inst_grid_name_imag
             inst_grid_real = readfits(inst_grid_name_real)
             inst_grid_imag = readfits(inst_grid_name_imag)
             inst_grid = dcomplex(inst_grid_real, inst_grid_imag)
             grid_input_structure = create_struct('inst_grid', inst_grid, grid_input_structure)
          endif
         
          grid_result = atm_and_inst_otf_grid(grid_input_structure, $
                                              corr = atm_corr, config_file=config_file, _Extra = extra)
          if save_otf EQ 1 THEN BEGIN
              print, 'FIND_STF: Saving OTF Grid Files'
              print, '   ', gridposFile
              print, '   ', otfgridFile_real
              print, '   ', otfgridFile_imag
              fits_write, otfgridFile_real, real_part(grid_result.atm_and_instr_grid)
              fits_write, otfgridFile_imag, imaginary(grid_result.atm_and_instr_grid)
              fits_write, gridposFile, grid_result.grid_fine
          endif
      
          ex = { atm_and_instr_grid:  ptr_new(grid_result.atm_and_instr_grid), $
                 grid:                ptr_new(grid_result.grid_fine) }
          undefine, atm_and_instr_grid
          undefine, grid
      endif else begin
          print, 'FIND_STF: Reading in OTF Grid Files'
          print, '   ', gridposFile
          print, '   ', otfgridFile_real
          print, '   ', otfgridFile_imag
          ;;;should read otfgridFile_real file before using the variable?
          grid = readfits(gridposFile)
          if not keyword_set(varPSF) then atm_and_instr_grid = complex(readfits(otfgridFile_real),readfits(otfgridFile_imag))
          ex = { atm_and_instr_grid:  ptr_new(atm_and_instr_grid), $
                 grid:                ptr_new(grid) }
          undefine, atm_and_instr_grid
          if not keyword_set(nearestneighbor) then undefine, grid 
       endelse
   
       if (N_ELEMENTS(extra) NE 0) then begin
          extra = CREATE_STRUCT(extra, ex)
       endif else begin
          extra = ex
       endelse
    endif

    if keyword_set(center_by_fitting) then begin
       extra = CREATE_STRUCT('center_by_fitting', center_by_fitting, extra)
    endif

    ;----------
    ; Load Input Starlist
    ;----------

    inStars = find_stf_lis_input_list(input)

    ;------------------------------
    ;
    ; Calculate noise and size of image
    ;
    ;------------------------------
    ; Default noise calc does median filtering... this shouldn't 
    ; be done on poorly sampled data (i.e. seeing limited data).
    if (scale GT 0.1) then nosub = 1 else nosub = 0

    print, format='(%"FIND_STF: Image noise statistics (nosub = %d)")', nosub
    gauss_noise_std, image, mode, std, patch=3, nterms_fit=6, nosub=nosub

    ;------------------------------
    ;
    ; Setup Starfinder paramters
    ;
    ;------------------------------
    if not keyword_set(quick) then begin
        psfIters = 3            ; Only used if makePsf = 1
        threshold = [3.,3.,3.]
        correl_mag = 2
        niter = 2
        n_fwhm_back = backBoxFWHM
    endif else begin
        print, "FIND_STF: Quick starfinder only."
        psfIters = 1            ; Only used if makePsf = 1
        threshold = [30.,30.]
        correl_mag = 2
        niter = 1
        n_fwhm_back = backBoxFWHM
     endelse

    x16c = inStars[0].x
    y16c = instars[0].y
    m16c = inStars[0].mag
    year = inStars[0].year
    
    ;------------------------------
    ;
    ; Get PSF
    ;
    ;------------------------------
    if (makePsf EQ 0) then begin
        ; Get background image. It should already exist if we
        ; also have a PSF.
        fits_read, bkgFile, background

        ; Use pre-existing PSF. Still need to repair saturated
        ; sources though.
        if not keyword_set(varPSF) then begin
            fits_read, psfFile, psf_2D, hdr_psf
            psf_fwhm = fwhm(psf_2D, MAG = 3, /CUBIC)
        endif else psf_fwhm = fwhm(varPSF[*,*,0], MAG = 3, /CUBIC)
        back_box = n_fwhm_back * psf_fwhm

        if aoopt eq 1 then begin
            if keyword_set(nearestneighbor) then begin
               if not keyword_set(partition) then begin
                  partition = nn_partition(grid[0,*], grid[1,*], siz[0])
               endif
               undefine, grid
            endif else begin
               partition = sf_partition(part_size, siz[0], over=grid_over)
            endelse
            
            if keyword_set(varPSF) then begin
               psf = ptr_new(varPSF)
            endif else begin
               psf = ptr_new(apply_otf_ratio(psf_2D, *(extra.atm_and_instr_grid)))
            endelse
            
            if tag_exist(extra,'debug') eq 1b then begin
               if extra.debug eq 1 then begin
                  writefits, 'psf_grid' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', *psf
               endif
            endif
            
            extra = CREATE_STRUCT('SV_PAR', partition, extra)
            if keyword_set(nearestneighbor) then begin
               extra = CREATE_STRUCT('nearestneighbor', nearestneighbor, extra)
            endif
            extra = CREATE_STRUCT(partition, extra) 
        endif

        ;; REPAIR saturated sources
        ;; Most of this code is hacked together from pieces of 
        ;; psf_extract.
        if not isSpeckle and exist(maxFile) then begin
           stf_psf_positions, image, x16c, y16c, year, scale, theta, $
                              xPsf, yPsf, xSec, ySec, $
                              filter=filter, /saturated, $
                              maxThresh=maxThreshold, refName=cooStar, $
                              starlist=starlist, psfSize=psf_size, $
                              box_psf=psfSearchBox, box_sec=secSearchBox
       
           ;; Find only saturated sources
           ;; Modified by S. Hornstein to find saturated sources where the
           ;; core has rolled over
           n_satur=0              
           peaks = image[xPsf, yPsf]

           for i=0L, n_elements(xPsf)-1 do $
              peaks[i]=max(image[xPsf[i]-2:xPsf[i]+2,yPsf[i]-2:yPsf[i]+2])

           w = where(peaks lt maxThreshold, count, $
                     complement=v, ncomplement=n_satur)

           if aoopt ne 1 then psf = ptr_new(psf_2D)

           if (n_satur GT 0) then begin
              x_sat = xPsf[v]  &  y_sat = yPsf[v]
          
              ; Fit and subtract secondary sources
              fitbox = round(2 * psf_fwhm) ; fitting box
              if aoopt eq 1 then begin
                 clean_image = psf_clean(image, psf, xSec, ySec, $
                                         fitbox, _EXTRA=extra)
              endif else begin
                 clean_image = psf_clean(image, psf, xSec, ySec, $
                                         fitbox, _EXTRA=extra)
              endelse
     
              repair_saturated, image, clean_image, background, psf_2D, $
                                psf_fwhm, x_sat, y_sat, maxThreshold, $
                                _EXTRA=extra
              if (makeRep NE 0) then fits_write, repFile, image, hdr
           endif
        endif else begin
          if not keyword_set(varPSF) then psf = ptr_new(psf_2D, /no_copy) else begin
             psf = ptr_new(varPSF)
             if n_elements(partition) eq 0 then partition = sf_partition(part_size, siz[0], over=grid_over)
             extra = CREATE_STRUCT('SV_PAR', partition, extra)
             extra = CREATE_STRUCT(partition, extra)
          endelse
       endelse

        xPsf = [x16c]
        yPsf = [y16c]

    endif else begin
        ;; First allow all saturated stars to be repaired by 
        ;; calling with the /saturated flag. Only bother with this
        ;; if it isn't speckle data (speckle has NO saturation)
        if not isSpeckle and exist(maxFile) then begin
            stf_psf_positions, image, x16c, y16c, year, scale, theta, $
                               xPsf, yPsf, xSec, ySec, $
                               filter=filter, /saturated, $
                               maxThresh=maxThreshold, refName=cooStar, $
                               starlist=starlist, psfSize=psf_size, $
                               box_psf=psfSearchBox, box_sec=secSearchBox

            ; Call the Starfinder PSF extraction code
            psf_extract, xPsf, yPsf, xSec, ySec, image, psf_size, $
                         psf_2D, psf_fwhm, background, $
                         N_FWHM_BACK=n_fwhm_back, iter=2, $
                         upper_level=maxThreshold, $
                         /rad_norm, _EXTRA=extra, N_SIGMA = N_SIGMA_1, MASKRAD = MASKRAD_1
                     
            if (makeRep NE 0) then fits_write, repFile, image, hdr
     
            back_box = n_fwhm_back * psf_fwhm
        endif

        ; Now re-extract throwing out some bad (saturated) PSF sources 
        ; such as IRS 13. They are probably extended or too crowded.
        stf_psf_positions, image, x16c, y16c, year, scale, theta, $
                           xPsf, yPsf, xSec, ySec, filter=filter, $
                           refName=cooStar, starlist=starlist, psfSize=psf_size, $
                           box_psf=psfSearchBox, box_sec=secSearchBox

        ; Call the Starfinder PSF extraction code
        psf_extract, xPsf, yPsf, xSec, ySec, image, psf_size, $
                     psf_2D, psf_fwhm, background, $
                     N_FWHM_BACK=n_fwhm_back, iter=2, $
                     upper_level=maxThreshold, $
                     /rad_norm, _EXTRA=extra, N_SIGMA = N_SIGMA_1, MASKRAD = MASKRAD_1
                 
        back_box = n_fwhm_back * psf_fwhm
        print, format='(%"FIND_STF: PSF_FWHM = %8.3f")', psf_fwhm
    
        ;; PSF fixes -- Cos-filter, Threshold, and Normalize
        if (fixPsf EQ 1) then begin
            if not keyword_set(quick) then begin
                psf_2D = fix_psf(psf_2D, psf_size, notrim=noPsfTrim, _Extra=extra)
            endif else begin
                psf_2D = fix_psf(psf_2D, psf_size, notrim=noPsfTrim, $
                                 innerHalo=2.0*psf_fwhm, outerHalo=3.0*psf_fwhm)
            endelse
        endif
    
        if tag_exist(extra,'debug') eq 1b then if extra.debug eq 1 then $
            writefits, 'fixed_psf_2D' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', psf_2D
         
        if (keyword_set(aoopt)) then begin
            if keyword_set(partition) then begin
                if keyword_set(nearestneighbor) then $
                    extra = CREATE_STRUCT('nearestneighbor', nearestneighbor, extra)
            endif else partition = sf_partition(part_size, siz[0], over=grid_over)
            psf = ptr_new(apply_otf_ratio(psf_2D, *(extra.atm_and_instr_grid)))  
            if tag_exist(extra,'debug') eq 1b then if extra.debug eq 1 then $
                writefits, 'psf_grid' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', *psf
            extra = CREATE_STRUCT('SV_PAR', partition, extra)
            extra = CREATE_STRUCT(partition, extra) 
        endif else psf = ptr_new(psf_2D)
    endelse
    
    ;--------------------
    ;   Clean up starlist and calculate fluxes.
    ;--------------------
    ; Fix bad magnitudes (they were wrong anyhow)
    idx = where(inStars.mag GT 25, cnt)
    if cnt GT 0 then inStars[idx].mag = 25.0

    ; Convert magnitudes to fluxes using coordinates for 
    ; the first star in the list and scaling from that. 
    ; Save the coordinates for later use.
    ; Convert all magnitudes to intensities
    f0 = max(image[x16c-4:x16c+4, y16c-4:y16c+4])
    inStars.i0 = (10^( -1 * (inStars.mag - m16c) / 2.5 )) * f0

    ; We need to trim out all the sources that are not
    ; within the size of this image. Lets also trim out all sources
    ; are just zeros in this image
    good = where(inStars.x GT 0+5 AND inStars.x LT siz[0]-5-1 AND $
                 inStars.y GT 0+5 AND inStars.y LT siz[1]-5-1, goodCnt)
    inStars=inStars[good]

    ;--------------------
    ; PERFORMANCE BOOST:
    ; Loop through and get rid of those with zeros
    ;--------------------
    boxSum = fltarr(n_elements(good))
    for i=0, n_elements(good)-1 do begin
        xl = round(inStars[i].x-2)
        xu = round(inStars[i].x+2)
        yl = round(inStars[i].y-2)
        yu = round(inStars[i].y+2)

        boxSum[i] = total(image[xl:xu,yl:yu])
     endfor

    good = where(boxSum NE 0, goodCnt)
    inStars=inStars[good]

    ; Count up the number of stars:
    n_max = goodCnt
    fmt = '(%"Candidates %d out of %d (16C @ %3d %3d with '
    fmt = fmt + 'mag=%4.1f flux=%8.1f)")'
    print, format=fmt, n_max, n_elements(inStars.x), x16c, y16c, m16c, f0

    x = inStars.x
    y = inStars.y
    f = inStars.i0

    ;;;------------------------------
    ;;;
    ;;; Loop through all the correlation values and run starfinder
    ;;;
    ;;;------------------------------
    for j=0, n_elements(corr)-1 do begin
        min_correlation = corr[j]
        
        starfinder, image, psf, $
                    background = background, $
                    BACK_BOX = back_box,  $
                    threshold, REL_THRESHOLD = 1, /PRE_SMOOTH, $
                    NOISE_STD = std, min_correlation, $
                    CORREL_MAG = correl_mag, $
                    DEBLEND = deblend, N_ITER = niter, $
                    x, y, f, sx, sy, sf, c, metric=metric, $
                    STARS = stars, FORCE = force,$
                     _EXTRA = extra, fix_pos = fix_pos

        if tag_exist(extra,'debug') eq 1b then if extra.debug eq 1 then begin
            writefits, 'stars' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', stars
            writefits, 'res' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', image-stars
        endif
    
        ; Trim fake sources
        if (trimfake NE 0) then begin
            print, 'FIND_STF: Trimming fake sources'
            find_stf_cut_fakes, psf, x, y, f, sx, sy, sf, c, stars, siz, _Extra=extra
        endif

        ; Repeat PSF extraction
        if (makePsf NE 0 and not keyword_set(quick)) then begin
            for p=0, psfIters-1 do begin
            
                print, '##############################'
                print, 'FIND_STF: THIS IS PSF RE-EXTRACTION STEP NO ', p+1, ' of ', psfIters
                print, '##############################'

                ; Recreate the original PSF (at the full size)
                psf_2D = superpose_stars(image - background, $
                                         xPsf, yPsf, psf_size, psf_fwhm, $
                                         /MAX_NORM, INTERP='I', $
                                         saturation=maxThreshold, _EXTRA=extra, $
                                         stars=stars, psf=psf, x_stars=x, y_stars=y, $
                                         fluxes = f, weighted = weighted)
                
                if tag_exist(extra,'debug') eq 1b then if extra.debug eq 1 then $
                    writefits, 'psf_2D_unclipped' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', psf_2D
                
                ;; if total(psf_2D) lt 0 then begin
                ;;     print, 'FIND_STF: WARNING setting negative PSF background to zero'
                ;;             ; this needs to be changed to really
                ;;             ; only subtract negative. Also
                ;;             ; retest.
                ;;             ; This was originally here to fix
                ;;             ; negative normalizations that occured
                ;;             ; when large numbers of halo pixels
                ;;             ; were negative. BUT perhaps we shoud
                ;;             ; ljust try adding a background
                ;;             ; constant?
                ;;             ; RETEST since secondary star fix was
                ;;             ; implemented. 
                ;;     psf_2D = cont_psf(psf_2D, n_sigma_1, maskrad_1)
                ;; endif
                ;;else psf_2D /= total(psf_2D)

                if tag_exist(extra,'debug') eq 1b then if extra.debug eq 1 then $
                    writefits, 'psf_2D_unclipped' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', psf_2D

                ; PSF fixes -- Extraction of contiguous PSF signal or Cos-filter, Threshold, and Normalize 
                ; Don't do cos and smoothing for AO datasets.
                ;; if back_clip then begin
                ;;    print, 'FIND_STF: Extracting contiguous PSF signal' 
                ;;    psf_2D = cont_psf(psf_2D, n_sigma_start - 0.1*float(p), maskrad)

                ;;    ; Attempt to smooth the PSF.
                ;;    ;;psf_2D = fix_psf(psf_2D, psf_size, nocos=noPsfCo, nosmooth=noPsfSmooth, _Extra=extra)
                ;;    ;; dist_circle, d_psf, 301   ; WHAT IS THIS HARDCODED 301? FIX!!!
                ;;    ;; m_psf = dblarr(301,301)
                ;;    ;; if field gt 0 then begin
                ;;    ;;    m_psf[where(d_psf lt 90.)] = 1.0d
                ;;    ;; endif else begin
                ;;    ;;    m_psf[where(d_psf lt 95.)] = 1.0d
                ;;    ;; endelse
                ;;    ;; m_psf = smooth(m_psf, 15) ; radial smoothing? Tapering to 0. Similar to fix_psf???
                ;;    ;; psf_2D = psf_2D * m_psf
                ;;    psf_2D /= total(psf_2D)
                ;;    print, 'total PSF flux: ' + string(total(psf_2D))
                ;;             ; Note this might not actually apply
                ;;             ; anymore because the GC clipping is
                ;;             ; already doing well???

                ;; endif else begin
                ;;     ;; psf_2D = fix_psf(psf_2D, psf_size, $
                ;;     ;;                  nocos=noPsfCo, nosmooth=noPsfSmooth, cliphalo=clipHalo, $
                ;;     ;;                  _Extra=extra)
                ;;     ;; PSF fixes -- Cos-filter, Threshold, and Normalize 
                ;;     ;; Don't do cos and smoothing for AO datasets.
                ;;     psf_2D = fix_psf(psf_2D, psf_size, nofix=noPsfFix)
                ;; endelse
                
                ;; PSF fixes -- Extraction of contiguous PSF signal or
                ;;              Cos-filter, Threshold, and Normalize
                ; Don't do cos and smoothing for AO datasets, only clipping.
                psf_2D = fix_psf(psf_2D, psf_size, nocos=noPsfCos, nosmooth=noPsfSmooth, $
                                 nohaloclip=noPsfClip, notrim=noPsfTrim, $
                                 clipsigma=n_sigma_start * (1.0 - (float(p) / psfIters)), $
                                 clipradius=maskrad)
            
                if tag_exist(extra,'debug') eq 1b then if extra.debug eq 1 then $
                    writefits, 'fixed_psf_2D' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', psf_2D
            
                ptr_free, psf
                if (keyword_set(aoopt)) then begin
                
                    psf = ptr_new(apply_otf_ratio(psf_2D, *(extra.atm_and_instr_grid)))
                
                    if tag_exist(extra,'debug') eq 1b then if extra.debug eq 1 then $
                        writefits, 'psf_grid' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', *psf
                
                endif else psf = ptr_new(psf_2D)

                starfinder, image, psf, $
                            background = background, $
                            BACK_BOX = back_box,  $
                            threshold, REL_THRESHOLD = 1, /PRE_SMOOTH, $
                            NOISE_STD = std, min_correlation, $
                            CORREL_MAG = correl_mag, $
                            DEBLEND = deblend, N_ITER = niter, $
                            x, y, f, sx, sy, sf, c, metric=metric, $
                            STARS = stars, force = force, _EXTRA = extra, fix_pos = fix_pos
            
                if tag_exist(extra,'debug') eq 1b then if extra.debug eq 1 then begin
                    writefits, 'stars' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', stars
                    writefits, 'res' + cgtimestamp(8,random_digits = 6, /valid) + '.fits', image-stars
                endif
            
                ;; Trim fake sources
                if (trimfake NE 0) then begin
                    print, 'FIND_STF: Trimming fake sources'
                    find_stf_cut_fakes, psf, x, y, f, sx, sy, sf, c, stars, siz, metric=metric, _Extra=extra
                endif

            endfor
         endif
        
        ;--------------------
        ; Do one more pass at starfinder (forced) to refine positions 
        ; and photometry if we have trimmed fakes.
        ;--------------------
        if (trimfake NE 0) then begin
           starfinder, image, psf, $
                       background = background, $
                       BACK_BOX = back_box,  $
                       threshold, REL_THRESHOLD = 1, /PRE_SMOOTH, $
                       NOISE_STD = std, min_correlation, $
                       CORREL_MAG = correl_mag, $
                       DEBLEND = deblend, N_ITER = niter, $
                       x, y, f, sx, sy, sf, c, metric=metric, $
                       STARS = stars, FORCE = 1, _EXTRA = extra, fix_pos = fix_pos
        endif
        
        out = [transpose(x), $
               transpose(y), $
               transpose(f), $
               transpose(sx), $
               transpose(sy), $
               transpose(sf), $
               transpose(c)]
    
        if (deblend EQ 0) then begin
            file = outRoot + '_' + strmid(strtrim(corr[j],1),0,3) 
        endif else begin
            file = outRoot + '_' + strmid(strtrim(corr[j],1),0,3)+'d'
        endelse

        ; Write output files
        if (makePsf NE 0) then fits_write, psfFile, psf_2D, hdr
        if (makePsf NE 0) then fits_write, bkgFile, background, hdr
        if (makeRes NE 0) then fits_write, resFile, image-stars-background, hdr
        if (makeStars NE 0) then fits_write, starsFile, stars, hdr
        if (aoopt NE 0) then fits_write, psfgridFile, *psf, hdr
        if ((aoopt NE 0) and keyword_set(makeGrid)) then fits_write, (inRoot + '_grid_pos.fits'), grid_result.grid_fine
    
        ;---------------
        ; Find 16C and make it the first star in the list
        ;---------------
        ; Find the stars closest to 16C with the function distance()
        ; then sort by flux and take the brightest, closest star.
        d = distance(x16c,y16c,x,y)
        w = where(d le 10.0)

        if n_elements(w) LT 0 then begin
            print, 'FIND_STF: coo star not detected in image!'
            print, 'Check coordinates in .coo file.'
        endif

        ; Sort based on flux
        fidx = reverse(sort(f[w]))
        w = w[fidx[0]] 

        ; Swap so that 16C is at the top.
        swap1 = out[*,0] & swap2 = out[*,w[0]]
        out[*,w[0]] = swap1 & out[*,0] = swap2

        ;----------
        ; Output files
        ;----------
        txtfile = file + '_stf.txt'
        OPENW, lun, txtfile, /GET_LUN
        PRINTF, lun, out  &  free_lun, lun
        OPENW, lun, file + '_metrics.txt', /GET_LUN
        PRINTF, lun, transpose(metric)  &  free_lun, lun

        txt2lis_counts, file+'_stf', sigFile, year, firstName=cooStar
        
        ; Trimframes algorithm. 
        if (keyword_set(trimframes)) then begin
            trimframes, file+'stf.lis'
        endif
    endfor

    ptr_free, psf

    if (keyword_set(aoopt)) then begin
        ptr_free, extra.atm_and_instr_grid
        ptr_free, extra.grid
    endif
end

