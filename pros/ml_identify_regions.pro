; Ranks highly correlated pixels, ensuring that none closer than
; width/2 are chosen
function ml_identify_regions, corr_vec, width

if width mod 2 eq 0 then message, "Width must be odd"

corr_vec_cp = corr_vec

corrpix = []

npix = n_elements(corr_vec_cp) 
npix_remain = npix

while npix_remain gt 0 do begin
    
    ; Find the pixel with the maximum magnitude of correlation
    max_corr = max(abs(corr_vec_cp), max_corr_idx, /nan)

    ; Save the index of the highly correlated pixels
    corrpix = [corrpix, max_corr_idx]

    ; NaN the indices in a window of width pixels, surrounding the
    ; highly correlated pixel
    remove_idx = lindgen(width)-(width/2) + max_corr_idx

    ; make sure indices are in range of vector
    masked_idx = where(remove_idx lt 0 or remove_idx ge npix, nmask)
    if nmask ne 0 then remove, masked_idx, remove_idx
    ; mask the window
    corr_vec_cp[remove_idx] = !values.d_nan

    ; Report the remaining pixels
    npix_remain = total(finite(corr_vec_cp))
    
endwhile

return, corrpix

end
