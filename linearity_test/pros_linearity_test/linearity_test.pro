;+
; NAME:
;
;
;
; PURPOSE:
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

function broaden_spectrum, spectrum, vsini, WL_GRID=wl_grid, EPSILON=epsilon, OVERSAMP=oversamp

if n_elements(epsilon) eq 0 then epsilon = 0.6
if n_elements(oversamp) eq 0 then oversamp = 11L

; FIXME: Make sure spectrum is on evenly spaced logwl grid

; Oversample spectrum
; FIXME: Maybe don't linearly interpolate the logwl grid (didn't stop
; me before...)
x_over = (dindgen(n_elements(wl_grid)*oversamp)-(oversamp/2L))/oversamp

oversampled_spectrum = interpol(spectrum, dindgen(n_elements(wl_grid)), x_over)

; Set up delta v
dlogwl = 6d-6
bigc = !const.c*1d-3  ;c in km/s
deltav = ((10d0^(dlogwl/oversamp))-1d0)*bigc

; Create rotation kernel
lsf_rot_vel = lsf_rotate(deltav, vsini, velgrid=velgrid, epsilon=epsilon)

; Convert velocity-space x-axis to logwl-space
velgrid_kernel_x = alog10(1d0 + (velgrid/bigc))

; Make a similar vector which matches the logwl-grid
nlsfrot = n_elements(velgrid)+6L ;pad to make sure we get the whole kernel
lsf_rot_kernel_x = (dindgen(nlsfrot)-(nlsfrot/2L)) * dlogwl/oversamp
; Spline lsf_rot onto the logwl_grid vector
lsf_rot = interpol(lsf_rot_vel, velgrid_kernel_x, lsf_rot_kernel_x)

; Make sure interpolation doesn't introduce negatives
rot_neg_idx = where(lsf_rot lt 0d0, ncnt)
if ncnt gt 0 then lsf_rot[rot_neg_idx]=0d0

; Normalize
norm_term = int_tabulated(lsf_rot_kernel_x, lsf_rot, /double)
;lsf_rot = lsf_rot/total(lsf_rot, /double)
lsf_rot = lsf_rot/norm_term

; Convolve with the spectrum in order to broaden
rot_spec = convol(oversampled_spectrum, lsf_rot, /edge_wrap)

; Downsample
broadened_spectrum = downsample_tophat(rot_spec, oversamp)

return, broadened_spectrum


end


pro linearity_test

; Readin spectra

; Find one with Vsini_CKS < 1

; Loop and broaden a bunch of times

; Find correlated pixels and plot slope vs. vsini

npix = 101L
wl_grid = dindgen(npix)*6d-6 + 4d0
x = dindgen(npix)-50

fakespec = replicate(1d0, npix)

for linenum = 0, 20 do begin

    depth = randomu(seed, 1)*0.5d0
    centroid = randomu(seed,1)*(npix-1) - (npix/2)
    
    fakespec = fakespec - gaussian(x, [depth, centroid, 1d0])

endfor

nvsini = 101
vsini_vec = dindgen(nvsini)

slope_result = dblarr(npix, nvsini)

foreach vsini, vsini_vec, idx do begin

    synth_spec = broaden_spectrum(fakespec, vsini, wl_grid=wl_grid)

    if idx eq 0 then synth_spec = fakespec

    deriv_spec = ml_savgol_spectra(wl_grid, synth_spec, width=5)

    slope_result[*,idx] = deriv_spec

endforeach

corr_vec = dblarr(npix)

foreach pixel, lindgen(npix), pixidx do begin

    correlation = correlate(reform(slope_result[pixel,*]), vsini_vec, /double)

    corr_vec[pixel] = correlation

endforeach

corrsort = sort(abs(corr_vec))

; Going from most correlated, plot

for i=0,100 do begin

    pixnum=i
    ;pixnum = corrsort[i]

    slope_at_pixel = reform(slope_result[pixnum,*])

    fitcoeff = robust_poly_fit(vsini_vec, slope_at_pixel, 5, polyfit, /double)

    
    plot, vsini_vec, slope_at_pixel, ps=8, xs=2, ys=2, tit="Slope at pixel= "+strtrim(i,2)

    oplot, vsini_vec, polyfit, /thick

    dxstop
    dummy=0

endfor


stop



end
