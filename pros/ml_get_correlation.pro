; Calculate correlations between columns of input spec_arr and elements
; of yparam
function ml_get_correlation, spec_arr, yparam

; Initialize
npix = n_elements(spec_arr[*,0]) 
corr_vec = dblarr(npix)

; Loop through and get correlation at each pixel
for i = 0, npix-1 do begin

    input_vec = spec_arr[i,*]
    
    corr = correlate(input_vec, yparam, /double)
    
    corr_vec[i] = corr

endfor

return, corr_vec

end
