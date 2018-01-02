; Make filelist for all the CKS-type dwarfs in APOGEE
pro apgbs_allfiles

; Select all stars with GKd or Fd classes,
; Teff between 4640 and 7140,
; FeH between -1.15 and 0.47

allapg = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allStar-l31c.2.fits',1)

; match class
class_match = stregex(allapg.aspcap_class, 'G*[FK]d_.*', /boolean) 

; Take only apo25m
apo25m_match = strmatch(strtrim(allapg.telescope,2), 'apo25m')

; Take only stars without STAR_BAD aspcapflag
flag_match = stregex(allapg.aspcapflags, '.*STAR_BAD.*', /boolean)
noflag_match = ~flag_match

; Take only stars with the right temperature
teff = allapg.teff
teff_match = (teff le 7140) and (teff ge 4640)

; Take only stars with the right metallicity
feh = allapg.m_h
feh_match = (feh ge -1.15) and (feh le 0.5)

; Make a SNR cut of >100
snr = allapg.snr
snr_match = snr ge 100


; MASTER MASK
match_vec = class_match and apo25m_match and noflag_match and teff_match and $
  feh_match and snr_match

match_idx = where(match_vec eq 1, nmatch)

; MASTER STRUCT
apg = allapg[match_idx]


;;;;; MAKE DOWNLOAD BATCH FILE
; Download the APG spectra - make list of filenames
loc = strtrim(apg.location_id,2)
filename = strtrim(apg.file,2)

idx = strtrim(match_idx,2)

full_list = idx + '|' + loc + '/' + filename

; Write the filenames to a text file
openw, lun, '/home/stgilhool/APGBS/spectra/apgbs_filelist_all_withidx.txt', /get_lun

foreach line, full_list, idx do begin

    printf, lun, line

endforeach

free_lun, lun
; Get the spectra with the following command (without spider)
;wget --spider -nv -r -nH --cut-dirs=8 \
;    -i apgbs_filelist_all.txt \
;    -B https://data.sdss.org/sas/dr14/apogee/spectro/redux/r8/stars/apo25m/

; Also write a structure that is the allStar.fits file, but with just
; the entries we want
mwrfits, apg, '/home/stgilhool/APGBS/spectra/allStar_selection.fits', /create


stop

end


