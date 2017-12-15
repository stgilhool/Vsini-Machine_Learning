; Quick one-off to add the ALLSTAR_IDX field to the file
; allStar_selection.fits.  The ALLSTAR_IDX is the index in the master
; allStar file, from which the corresponding entry of allStar_slice.fits was taken
pro add_allstaridx

; Readin
; this one hold allstar_idx, but also spectra, so it's slow to read
alldata_file = '/home/stgilhool/APGBS/spectra/apgbs_datacube.fits'
ad_full = mrdfits(alldata_file, 1)

; this is our sample's subset of the allStar file
allstar_file = '/home/stgilhool/APGBS/spectra/allStar_slice.fits'
astar = mrdfits(allstar_file,1)

n_data = n_elements(ad_full)  
n_astar = n_elements(astar) 

if n_data eq n_astar then print, "N_elements agree" $
  else print, "N_elements disagree"


for i = 0, n_data-1 do begin

    ; grab the index
    allstar_idx = ad_full[i].allstar_idx

    ; grab the allstar info structure
    astruct = astar[i]

    ; create a new structure
    newstruct = create_struct('ALLSTAR_IDX', allstar_idx, astruct)

    ; save
    if i eq 0 then outstruct = replicate(newstruct, n_data)
    outstruct[i] = newstruct

endfor

; TEST
; help, outstruct

; randomidx = fix(randomu(seed, 20)*(n_data-1))

; outtest = outstruct[randomidx]
; intest = ad_full[randomidx]

; printcol, intest.allstar_idx, outtest.allstar_idx, intest.teff_apg, outtest.teff

; stop

; WRITE
mwrfits, outstruct, '/home/stgilhool/APGBS/spectra/allStar_selection.fits', /create


end
