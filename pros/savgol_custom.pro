;+
; NAME: SAVGOL_CUSTOM
;
;
;
; PURPOSE: To perform Savitsky-Golay filtering, with the option of
; returning errors on SG parameters
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE: savgol_custom(wl_grid, spec_arr, err_arr,
;                                 WIDTH=width, ORDER=order,
;                                 DEGREE=degree, 
;                                 SAVGOL_ERROR=savgol_error)
;
;
;
; INPUTS:
;    wl_grid - the evenly spaced wavelength vector 
;    spec_arr - the npix X nspectra data vector
;    err_arr - the npix X nspectra flux error vector
;
; OPTIONAL INPUTS:
;    
;
;
; KEYWORD PARAMETERS:
;    WIDTH - The width of the savgol filter window. Default = 11
;    ORDER - The desired order of returned data. Default = 1 (first
;            derivative)
;    DEGREE - The degree of the polynomial fit performed by savgol
;             filter. Default = 4
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;    SAVGOL_ERROR - The 1-sigma uncertainties on the returned data
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
;    WL_GRID is assumed to be evenly spaced. The function will permit
;    non-evenly spaced wl_grid, as long as the delta_wl doesn't change
;    by more than 0.1 * mean(delta_wl)
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;    To perform ordinary savgol smoothing:
;       smoothed_spectrum = savgol_custom(wl_grid, spec_arr, err_arr,
;       ORDER=0)
;
;    To get smoothed 1st derivative, 5-pixel window and error
;    estimates:
;        smooth_deriv = savgol_custom(wl_grid, spec_arr, err_arr,
;        WIDTH=5, SAVGOL_ERROR=savgol_error)
;
; MODIFICATION HISTORY:
;
;       Mon Jan 15 11:51:37 2018,
;       <stgilhool@iroquois.physics.upenn.edu>
;
;		
;
;-

function savgol_custom, wl_grid, spec_arr, err_arr, WIDTH=width, ORDER=order, DEGREE=degree, SAVGOL_ERROR=savgol_error

savgol_error = 0

; Keywords
if n_elements(width) eq 0 then width = 11
if n_elements(order) eq 0 then order = 1
if n_elements(degree) eq 0 then degree = 4

; prep wl_grid (automatically take average delta x)
dx2 = deriv(dindgen(n_elements(wl_grid)), wl_grid)
deltax = mean(dx2)

; if delta x changes a lot in the supplied wl_grid, halt execution
if max(dx2)-min(dx2) gt 0.1*deltax then message,"wl grid is too non-linear" else $
  if max(dx2)-min(dx2) ne 0 then begin
    print, "WARNING: wl_grid is not evenly spaced! Result will be biased!"
    wait, 1
endif

; nspectra and npix
spec_arr_size = size(spec_arr, /dim)

npix = spec_arr_size[0]

if n_elements(spec_arr_size) eq 1 then nspectra = 1 else $
  if n_elements(spec_arr_size) eq 2 then nspectra = spec_arr_size[1] else $
  message, "spec_arr must be 1 or 2 dimensional"
    
; Simulate the effect of the filter

; Create design matrix
z = dindgen(width) - (width/2)

;bigm = replicate(1d0, degree+1, width)
bigm = dblarr(degree, width)

; populate columns of design matrix
for degnum = 1, degree do begin

    bigm[(degnum-1),*] = z^degnum

endfor

; Create array to hold results
coeff_results = dblarr(npix, nspectra, degree+1)
error_results = dblarr(npix, nspectra, degree+1)

; Loop through spectra
for snum = 0, nspectra-1 do begin

    tempspec = spec_arr[*,snum]
    temperr = err_arr[*,snum]

    ; Pad the ends spec
    start_pad = reverse(tempspec[0:(width/2)-1])
    end_pad = reverse(tempspec[-1*(width/2):-1])
    tempspec = [start_pad, tempspec, end_pad]
    
    ; Pad the ends of err
    err_start_pad = reverse(temperr[0:(width/2)-1])
    err_end_pad = reverse(temperr[-1*(width/2):-1])
    temperr = [err_start_pad, temperr, err_end_pad]
    

    ; Step through all the pixels
    for pixnum = 0, npix-1 do begin

        ; make window index vector for tempspec (so no longer center window)
        window_idx = lindgen(width) + pixnum
        ; Take data in the window and find the coefficients
        window_data = tempspec[window_idx]
        window_err = temperr[window_idx]

        coeffs = regress(bigm, window_data, measure_errors=window_err, $
                         sigma=param_errors, const=rconst)
        
        coeffs = reform(coeffs)
        
        ; If regression fails, 0 the returned coefficients, and set
        ; errors high
        if total(finite([rconst,coeffs])) ne degree+1 then begin
            rconst = 0d0
            coeffs = replicate(0d0, degree)
            param_errors = replicate(1d8, degree)
        endif
        
        ; Store answers
        coeff_results[pixnum,snum,*] = [rconst,coeffs]
        error_results[pixnum,snum,*] = [0,param_errors]
    endfor
endfor

; Return smoothed spectra of desired order
normalization_term = factorial(order)/(deltax^order)

smooth_spec = coeff_results[*,*,order]*normalization_term
savgol_error = error_results[*,*,order]*normalization_term

return, smooth_spec
        
end
   
