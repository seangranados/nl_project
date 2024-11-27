pro starplant, d, int_num, rootdir, psfFile, epoch, dr_dir, imsize, filt, stf_lp_lis, input

    ; Plants one star with a given magnitude at a random placement from S0-2 given a distance
    ; Needs to be run in combo directory
    ;
    ; Parameters 
    ; -----------
    ; mag : int
    ;      magnitude of planted source
    ; d : int
    ;      distance from S0-2
    ; int_num : int
    ;      simulation run number for naming purposes (i.e. files name from MAG_DIST_{INT_NUM})
    ; star : str
    ;      name of star
    ; rootdir : str
    ;      location of directories above parent directory
    ; psfFile : str
    ;      PSF file to use for starplanting
    ; s02_x : int
    ;      x position of S0-2 in image (or star want to plant around)
    ; s02_y : int
    ;      y position of S0-2 in image (or star want to plant around)
    ; epoch : str
    ;      epoch planting in 
    ; dr_dir : str
    ;      location of image files or data release to grab images from
    ; counts : float
    ;      counts of desired magnitude star
    ; mag_ref : float
    ;      magnitude that corresponds to the number of counts (defined for each epoch)
    ; imsize : int
    ;      size of starplanting image 
    ; filt : str
    ;      filter of image

    seed = !NULL
	psf = readfits(psfFile)
	
    ; Count the number of stars to deal with
    numStars= file_lines(input)

    ; If something unexpected happens get out
    if (numStars lt 1) then begin
        print, "****** 0 stars in the input list. SOMETHING WRONG. ******"
        retall
     endif

    ; Create variables to store stuff
    inMag = fltarr(numStars)
    inX = fltarr(numStars)
    inY = fltarr(numStars)
    inF = fltarr(numStars)

    ;Read in the File
    openr, u, input, /get_lun

    ; Loop through the file and read in all the info
    img_gen = image_model(0, 0, 0, imsize, imsize, psf)
	for i=0, numStars-1 do begin
        temp = ''
        readf, u, temp
        line = strsplit(temp, /EXTRACT)

        ; These are the same in either case
        inMag[i] = float(line[1])
        inX[i] = float(line[2])
        inY[i] = float(line[3])
		inF[i] = float(line[4])
		
		img_gen = img_gen + image_model(inX[i], inY[i], inF[i], imsize, imsize, psf)
		print, line[0], ' ', line[2], ' ', line[3]
		
		
	endfor
		

    ; image to plant in from dr1
    im_file = rootdir+'/'+epoch+'/'+epoch+'_'+stf_lp_lis+'_model.fits'
    print, im_file
    img_prev = readfits(im_file, header)
    img = img_prev + img_gen

    ; saves new image with planted star in appropriate mag_dis_run directory
    savefile = rootdir+'/'+epoch+'/'+epoch+'_'+stf_lp_lis+'_model.fits'
    print, savefile
    savefile = STRCOMPRESS(savefile, /REMOVE_ALL)
    writefits, savefile, img, header

end