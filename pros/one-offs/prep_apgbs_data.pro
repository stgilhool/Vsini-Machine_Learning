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
;       Mon Oct 9 09:18:47 2017,
;       <stgilhool@iroquois.physics.upenn.edu>
;
;		
;
;-

pro prep_apgbs_data, FLUX=flux

if n_elements(flux) eq 0 then flux = 0

nx = 8575L


;;; Overlaps file
olstr = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/CKS_APOGEE_overlaps.fits',1)

;;; CKS data
cksinfo = mrdfits('/home/stgilhool/CKS/table/cks_table.fits',1)
cksinfo = cksinfo[olstr.tab_idx]

if strtrim(cksinfo[100].id_kic,2) ne strtrim(olstr[100].id_kic,2) then message, "CKS mismatch"

;;; allStar data
apginfo = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allStar-l31c.2.fits',1)
apginfo = apginfo[olstr.apg_idx]

if strtrim(apginfo[101].apogee_id,2) ne strtrim(olstr[101].id_2m,2) then message, "APG mismatch"

;;; APOGEE files
root_dir = '/home/stgilhool/APGBS/'
data_path = root_dir + 'spectra/overlaps/'

listfile = data_path + 'apgbs_filelist.txt'


; Readin file names
readcol, listfile, datafile, format="A", count=nfiles_tot

if nfiles_tot ne n_elements(apginfo) then message, "Wrong number of files?"

filenames = data_path + datafile

nfiles=nfiles_tot

; Read in headers

header = ptrarr(nfiles, /allocate_heap)

spectra = dblarr(nx,nfiles)
errors = dblarr(nx,nfiles)
mask = lonarr(nx,nfiles)

x = dindgen(nx)

; Mask flags
flag_bit = [0,1,2,3,4,5,6,12,14]
flag_decimal = long(total(2L^flag_bit))


;logwl_grid = dblarr(nx, nfiles)
;log_shift = dblarr(nfiles)

;for i = 0, nfiles-1 do begin
outstruct = []

foreach file, filenames, idx do begin

    ; Open fits file
    if file_test(file) then begin
        fits_open, file, fcb
        
        ; Read in header for a file
        fits_read, fcb, null, head, /header_only, exten_no=0
        
        ; save header
        *(header[idx]) = head
        
        ; Read in 1d spectrum array for a file
        fits_read, fcb, specall, exten_no=1
        fits_read, fcb, errall, exten_no=2
        fits_read, fcb, maskall, exten_no=3
        fits_close, fcb
        
        ; Save only the first (combined) spectrum
        spec = double(specall[*,0])
        err = double(errall[*,0])
        msk = maskall[*,0]
        
        mask_i = flag_decimal and msk
        
        ; normalize it
        if flux eq 0 then begin
            ; Remove continuum
            bmin = fxpar(head, 'BOVERMIN')
            bmax = fxpar(head, 'BOVERMAX')
            gmin = fxpar(head, 'GOVERMIN')
            gmax = fxpar(head, 'GOVERMAX')
            rmin = fxpar(head, 'ROVERMIN')
            rmax = fxpar(head, 'ROVERMAX')
            
            nblue  = bmax - bmin+1
            ngreen = gmax - gmin+1
            nred   = rmax - rmin+1
            
            bnorm = continuum_fit2(dindgen(nblue), spec[bmin:bmax], pix_mask=mask_i[bmin:bmax])
            gnorm = continuum_fit2(dindgen(ngreen), spec[gmin:gmax], pix_mask=mask_i[gmin:gmax])
            rnorm = continuum_fit2(dindgen(nred), spec[rmin:rmax], pix_mask=mask_i[rmin:rmax])
            
            ; Normalize and scale errors for each chip
            spec[bmin:bmax] = spec[bmin:bmax]/bnorm
            spec[gmin:gmax] = spec[gmin:gmax]/gnorm
            spec[rmin:rmax] = spec[rmin:rmax]/rnorm
            
            err[bmin:bmax] = err[bmin:bmax]/bnorm
            err[gmin:gmax] = err[gmin:gmax]/gnorm
            err[rmin:rmax] = err[rmin:rmax]/rnorm
            
            ; Zero the regions in between
            spec[0:bmin-1] = 0d0
            spec[bmax+1:gmin-1] = 0d0
            spec[gmax+1:rmin-1] = 0d0
            spec[rmax+1:-1] = 0d0
            
            err[0:bmin-1] = 0d0
            err[bmax+1:gmin-1] = 0d0
            err[gmax+1:rmin-1] = 0d0
            err[rmax+1:-1] = 0d0
            
            
        endif 
        
        ; Get vhelio, and wl coeffs for a file
        ; (vhelio is unnecessary, and wl coeffs are same for all)
        ;vhelio = fxpar(head, 'VHELIO', datatype=0.0d)
        wl0    = fxpar(head, 'CRVAL1', datatype=0.0d)
        wl1    = fxpar(head, 'CDELT1', datatype=0.0d)
        
        spectra[*,idx] = spec
        errors[*,idx] = err
        mask[*,idx] = msk
        
    endif else begin
        
        spec = replicate(!values.d_nan, nx)
        err = replicate(!values.d_nan, nx)
        msk = replicate(0, nx)

    endelse
    
    outstr1 = {OVERLAP_IDX:idx, $
               APG_IDX:olstr[idx].apg_idx, $
               CKS_IDX:olstr[idx].tab_idx, $
               TEFF_APG:apginfo[idx].teff, $
               TEFF_CKS:cksinfo[idx].iso_steff, $
               FEH_APG:apginfo[idx].m_h, $
               FEH_CKS:cksinfo[idx].iso_smet, $
               LOGG_APG:apginfo[idx].logg, $
               LOGG_CKS:cksinfo[idx].iso_slogg, $
               VSINI_APG:apginfo[idx].vsini, $
               VSINI_CKS:cksinfo[idx].cks_svsini, $
               SPEC:spec, $
               ERR:err, $
               MASK:msk, $
               MASK_BIT:flag_decimal}
    
    outstruct = [outstruct, outstr1]
    
endforeach

; Make logarithmic wl scale and logarithmic wl shift
logwl_grid = wl0 + (x * wl1)

outfile = '/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits'

; Write the wl_grid to the HDU0
mwrfits, logwl_grid, outfile, /create
; Write data structure to HDU1
mwrfits, outstruct, outfile

stop

end
