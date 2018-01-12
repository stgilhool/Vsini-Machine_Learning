; Generate a Savitsky-Golay filter and apply it to spectra data
function ml_savgol_spectra, wl_grid, spec_arr, WIDTH=width, ORDER=order, DEGREE=degree

; Keywords
if n_elements(width) eq 0 then width = 11
if n_elements(order) eq 0 then order = 1
if n_elements(degree) eq 0 then degree = 4

; prep wl_grid (automatically take average delta x)
dx2 = deriv(dindgen(n_elements(wl_grid)), wl_grid)
deltax = mean(dx2)

; if delta x changes a lot in the supplied wl_grid, halt execution
if max(dx2)-min(dx2) gt 0.1*deltax then message,"wl grid is too non-linear"

; nspectra
spec_arr_size = size(spec_arr, /dim)

if n_elements(spec_arr_size) eq 1 then nspectra = 1 else $
  nspectra = spec_arr_size[1]

; Apply savgol filter
nleft = width/2
nright = width/2

; make filter
filter = savgol(nleft, nright, order, degree, /double)
sgfilter = (factorial(order)/(deltax^order)) * filter

; apply filter
smooth_spectra = spec_arr
for i = 0, nspectra-1 do begin

    tempspec = spec_arr[*,i]

    ; interpolate over masked (NaNed) pixels
    interp_spec = interpol(tempspec, wl_grid, wl_grid, /nan)

    tempspec = interp_spec

    ; apply the filter through convolution
    smoothspec = convol(tempspec, sgfilter, /edge_wrap)

    ; Replace the spectrum in smooth_spectra with the smoothed version
    smooth_spectra[*,i] = smoothspec

endfor

return, smooth_spectra

end

