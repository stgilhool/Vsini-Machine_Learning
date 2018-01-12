;+
; NAME: SLOPE_ERROR_TEST
;
;
;
; PURPOSE: To test a couple of ways of determining slope error
;
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
if max(dx2)-min(dx2) gt 0.1*deltax then message,"wl grid is too non-linear"

; nspectra and npix
spec_arr_size = size(spec_arr, /dim)

npix = spec_arr_size[0]

if n_elements(spec_arr_size) eq 1 then nspectra = 1 else $
  if n_elements(spec_arr_size) eq 2 then nspectra = spec_arr_size[1] else $
  message, "spec_arr must be 1 or 2 dimensional"

; Simulate the effect of the filter

; Create design matrix
z = dindgen(width) - (width/2)

bigm = replicate(1d0, degree+1, width)

print, z, format='("Z= ", I)'

help, bigm
print, bigm

help, width
print, width/2, width, format='("Half-width= ", I, " for W= ", I)

stop

; populate columns of design matrix
for degnum = 1, degree do begin

    bigm[degnum,*] = z^degnum

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

        coeffs = regress(bigm[1:*,*], window_data, measure_errors=window_err, $
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
    
        

pro slope_error_test

; Read in a test spectrum
data = mrdfits('slope_err_test_spec.fits',1)
logwl_grid = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',0)

mask = data.mask
spec = data.spec
err = data.err

; Apply mask
bad_pix = where(mask ne 0, nbad)

err[305] = 1d0
;err[bad_pix] = 1d8

; Apply SAVGOL filter
slope_savgol = ml_savgol_spectra(logwl_grid, spec, width=5)

; Try manual SAVGOL approach to see if that can return errors
slope_custom = savgol_custom(logwl_grid, spec, err, width=5, savgol_error=savgol_error)

plot, slope_savgol[300:1000], /thick
oplot, slope_custom[300:1000], /thick, co=!red

slope = deriv(logwl_grid, spec)

printcol, slope_savgol[300:310], slope_custom[300:310], (spec[301:311]-spec[299:309])/(2*(6d-6)), slope[300:310], err[300:310], savgol_error[300:310]

; okay, so this appears to work. One cool thing is that for a single
; pixel with high error, the slope error is not bad. This makes sense
; since the slope is really calculated from the neighboring points.
; It's like an automatic interpolation with errors incorporated for
; free.

; I do have to set the error not too high, or else NaNs pop up. I
; think err=1 is actually perhaps high enough for all masked pixels,
; because we do know that only 0-1 are physical in a continuum
; normalized spectrum (modulo noise, so maybe 0-1.5).  So, assuming a
; decent continuum normalization, it's not to bad to set the flux
; error to only 1 or 2 or something, rather than 1d8. The errors for
; neighboring pixels still blow up, so I just have to make sure to set
; the flux error high enough that the errors blow up enough for me to
; flag them as bad slope pixels.  Then we can set the slope error to
; some high level.
stop
; Determine slope error from analytical thing



; Try Monte Carlo approach

end
