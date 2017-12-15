;+
; NAME:
;
;
;
; PURPOSE:
;  To read APOGEE spectra into a structure for easy handling
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;       Mon Oct 23 09:18:40 2017,
;       <stgilhool@iroquois.physics.upenn.edu>
;
;		
;
;-

pro readin_apgbs, NFILES=nfiles, FILES=files, FLUX=FLUX

if n_elements(nfiles) eq 0 and n_elements(files) eq 0 then begin
    nfiles=50
    filelist = lindgen(nfiles)
endif else if n_elements(nfiles) gt 0 and n_elements(files) gt 0 then begin
    print, "Ignoring NFILES, using FILES..."
    filelist = files
    nfiles = n_elements(filelist) 
endif else if n_elements(nfiles) gt 0 and n_elements(files) eq 0 then begin
    filelist = lindgen(nfiles)
endif else if n_elements(nfiles) eq 0 and n_elements(files) gt 0 then begin
    filelist = files
    nfiles = n_elements(filelist) 
endif


if n_elements(flux) eq 0 then flux = 0

bigc = 299792.458 ;km/s
nx = 8575L

root_dir = '/home/stgilhool/APGBS/'
data_path = root_dir + 'spectra/'

; Readin APG data
apginfo = mrdfits('/home/stgilhool/APGBS/spectra/allStar_selection.fits',1)

; Readin file names
listfile = data_path + 'apgbs_filelist_all_withidx.txt'

readcol, listfile, allstar_idx, datafile, format="L,A", delimiter='|', count=nfiles_tot

filenames = data_path + datafile

; Read in headers

;header = ptrarr(nfiles_tot, /allocate_heap)

;spectra = dblarr(nx,nfiles_tot)
;errors = dblarr(nx,nfiles_tot)
;mask = lonarr(nx,nfiles_tot)

x = dindgen(nx)

; Mask flags
flag_bit = [0,1,2,3,4,5,6,12,14]
flag_decimal = long(total(2L^flag_bit))


;outstruct = []

foreach file, filenames, idx do begin

    ; Open fits file
    if file_test(file) then begin
        fits_open, file, fcb
        
        ; Read in header for a file
        fits_read, fcb, null, head, /header_only, exten_no=0
        
        ; save header
        ;*(header[idx]) = head
        
        ; Read in 1d spectrum array for a file
        fits_read, fcb, specall, exten_no=1
        ;fits_read, fcb, errall, exten_no=2
        ;fits_read, fcb, maskall, exten_no=3
        fits_close, fcb
        
        ; Save only the first (combined) spectrum
        spec = double(specall[*,0])
        ;err = double(errall[*,0])
        ;msk = maskall[*,0]
        
        ;mask_i = flag_decimal and msk
        
        ; normalize it
        if flux eq 0 then begin
            ; Remove continuum
            bmin = fxpar(head, 'BOVERMIN')
            bmax = fxpar(head, 'BOVERMAX')
            gmin = fxpar(head, 'GOVERMIN')
            gmax = fxpar(head, 'GOVERMAX')
            rmin = fxpar(head, 'ROVERMIN')
            rmax = fxpar(head, 'ROVERMAX')
            
            ; Deal with weird problem
            if bmin gt bmax then begin
                bmin = fxpar(head, 'BMIN')
                bmax = fxpar(head, 'BMAX')
                bpix_flag = 1
            endif else bpix_flag = 0

            if gmin gt gmax then begin
                gmin = fxpar(head, 'GMIN')
                gmax = fxpar(head, 'GMAX')
                gpix_flag = 1
            endif else gpix_flag = 0
            
            if rmin gt rmax then begin
                rmin = fxpar(head, 'RMIN')
                rmax = fxpar(head, 'RMAX')
                rpix_flag = 1
            endif else rpix_flag = 0

            nblue  = bmax - bmin+1
            ngreen = gmax - gmin+1
            nred   = rmax - rmin+1
            
            

            bnorm = continuum_fit2(dindgen(nblue), spec[bmin:bmax]);, pix_mask=mask_i[bmin:bmax])
            gnorm = continuum_fit2(dindgen(ngreen), spec[gmin:gmax]);, pix_mask=mask_i[gmin:gmax])
            rnorm = continuum_fit2(dindgen(nred), spec[rmin:rmax]);, pix_mask=mask_i[rmin:rmax])
            
            ; Normalize and scale errors for each chip
            spec[bmin:bmax] = spec[bmin:bmax]/bnorm
            spec[gmin:gmax] = spec[gmin:gmax]/gnorm
            spec[rmin:rmax] = spec[rmin:rmax]/rnorm
            
            ;err[bmin:bmax] = err[bmin:bmax]/bnorm
            ;err[gmin:gmax] = err[gmin:gmax]/gnorm
            ;err[rmin:rmax] = err[rmin:rmax]/rnorm
            
            ; Zero the regions in between
            spec[0:bmin-1] = 0d0
            spec[bmax+1:gmin-1] = 0d0
            spec[gmax+1:rmin-1] = 0d0
            spec[rmax+1:-1] = 0d0
            
            ;err[0:bmin-1] = 0d0
            ;err[bmax+1:gmin-1] = 0d0
            ;err[gmax+1:rmin-1] = 0d0
            ;err[rmax+1:-1] = 0d0
            
            
        endif 
        
        ; Get wl solution
        wl0    = fxpar(head, 'CRVAL1', datatype=0.0d)
        wl1    = fxpar(head, 'CDELT1', datatype=0.0d)
        
        ;spectra[*,idx] = spec
        ;errors[*,idx] = err
        ;mask[*,idx] = msk
        
    endif else begin
        
        ;spec = replicate(!values.d_nan, nx)
        ;err = replicate(!values.d_nan, nx)
       ; msk = replicate(0, nx)

    endelse
    
    outstr1 = {APG_IDX:idx, $
               ALLSTAR_IDX:allstar_idx[idx], $
               TEFF_APG:apginfo[idx].teff, $
               FEH_APG:apginfo[idx].m_h, $
               LOGG_APG:apginfo[idx].logg, $
               VSINI_APG:apginfo[idx].vsini, $
               BPIX_FLAG:bpix_flag, $
               GPIX_FLAG:gpix_flag, $
               RPIX_FLAG:rpix_flag, $
               SPEC:spec};, $
               ;ERR:err, $
               ;MASK:msk, $
               ;MASK_BIT:flag_decimal}
    
    if idx eq 0 then outstruct = replicate(outstr1, nfiles_tot)

    outstruct[idx] = outstr1
    ;outstruct = [outstruct, outstr1]

    print, "Read File: "+strtrim(idx,2)+"/"+strtrim(nfiles_tot,2)
    
endforeach

; Make logarithmic wl scale and logarithmic wl shift
logwl_grid = wl0 + (x * wl1)

outfile = '/home/stgilhool/APGBS/spectra/apgbs_datacube.fits'

; Write the wl_grid to the HDU0
mwrfits, logwl_grid, outfile, /create
; Write data structure to HDU1
mwrfits, outstruct, outfile

stop

end
