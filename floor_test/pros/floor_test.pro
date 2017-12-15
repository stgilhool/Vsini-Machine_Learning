;+
; NAME: FLOOR_TEST
;
;
;
; PURPOSE: Run procedure to test if we can do better than
; resolution/instrumental-limited Vsini detection floor
;
;
;
; CATEGORY: Analysis
;
;
;
; CALLING SEQUENCE: floor_test, SPTYPE=sptype
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
;    SPTYPE - spectral type to run analsys on
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

pro floor_test, SPTYPE=sptype

;;; Check inputs
if n_elements(sptype) eq 0 then sptype_choice = 'default' $
  else sptype_choice = sptype

;;; First, grab some data
;Need spectra, parameters, and map between spectra and parameters
; need cks vsini and apogee spectral types

; Data for the overlaps
overlaps_data = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',1)
logwl_grid = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',0)
; Relevant CKS vectors
;cks_teff = overlaps_data.teff_cks
;cks_vsini = overlaps_data.vsini_cks
teff_ol_apg = overlaps_data.teff_apg
vsini_ol_cks = overlaps_data.vsini_cks

;;; Loop through temperature bins and print cks vsini range, and
;;; number of spectra

;t0_vec = maken(4650L, 6600L, (6600L-4650L)/100)
; bs = 100
; cks_teff_hist = histogram(cks_teff, reverse_indices=ri, min=4650L, binsize=bs)
; histx = findgen(n_elements(cks_teff_hist))*bs + 4650 


; plot, histx, cks_teff_hist, ps=10

t0 = 6000
t1 = 6170

sel_idx_ol = where(teff_ol_apg ge t0 and teff_ol_apg lt t1 and vsini_ol_cks le 25, $
                   nsel_ol)

odata = overlaps_data[sel_idx_ol]
teff_ol = odata.teff_apg
vsini_ol = odata.vsini_cks
spectra_ol = odata.spec

;;; Other data
;as = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allStar-l31c.2.fits',1)
; Read in APOGEE data cube
alldata_file = '/home/stgilhool/APGBS/spectra/apgbs_datacube.fits'
ad_full = mrdfits(alldata_file, 1)

; conditional statements defining temperature range 
c0 = 'allstar.teff_apg ge '+strtrim(t0,2)
c1 = 'allstar.teff_apg lt '+strtrim(t1,2)

nboot = n_bootstrap(c0, c1, allstar=ad_full, ret_idx=ri)

boot_astar_idx = ad_full[ri].allstar_idx
ol_astar_idx = odata.apg_idx

final_astar_idx = cgsetintersection(boot_astar_idx, ol_astar_idx, indices_a=bastar_idx, indices_b=oastar_idx)


; OKAY, you need to check the SNR cut? you made in making the data
; cube.  I think there are three stars in the OVERLAP set that aren't
; in the datacube.  The one I checked had SNR ~ 80 and no bad flags
; that I could see.  The others could possibly have bad flags, I can't
; remember if I checked that or not.

; Further pare down the overlap sample
odata = odata[oastar_idx]
teff_ol = teff_ol[oastar_idx]
vsini_ol = vsini_ol[oastar_idx]
spectra_ol = odata.spec

; Now, eliminate the overlaps from the bootstrap sample
final_datacube_idx = cgsetdifference(boot_astar_idx, final_astar_idx, positions=new_boot_idx)

boot_data = ad_full[ri[new_boot_idx]]

; test if that worked
help, ad_full[ri]
help, boot_data
if n_elements(ad_full[ri])-n_elements(boot_data) eq 48 then print, "N_ELEMENTS reduced by 48, as expected" else print, "Uh oh, overlap overlaps suspect!"
elim_test = cgsetintersection(boot_data.allstar_idx, final_astar_idx, success=elim_check)
if elim_check then print, "Uh oh, there are still overlaps!" else print, "Okay: no stars in the overlap sample present in the boot sample :)"

; Yay!


; sel_idx = where(ad_full.teff ge t0 and ad_full.teff lt t1, nsel)

; ad = ad_full[sel_idx]

; spectra = ad.spec

; ; Smooth all

;smooth_spec = ml_savgol_spectra(logwl_grid, spectra, width=5)


;apg_sdata = as[ad.allstar_idx]




; Partition
param_arr = [vsini_ol]

nsets = 2L

;divide points between training and cross-val
npts_arr = replicate(nsel_ol/nsets, nsets)
nmod = nsel_ol mod nsets
npts_arr[0] = npts_arr[0] + nmod

set_str = ml_partition_data(param_arr, npts_arr, /span)

train_idx = set_str.(0)
ntrain = n_elements(train_idx) 
cross_idx = set_str.(1)
ncross = n_elements(cross_idx) 

train_idx = lindgen(23)
cross_idx = lindgen(22) + 23

; Smooth that shit
smooth_spec_ol = ml_savgol_spectra(logwl_grid, spectra_ol, width=5)

; Correlations
corr_vec = ml_get_correlation(smooth_spec_ol, vsini_ol)

; choose regions
rwidth = 9
corrpix = ml_identify_regions(corr_vec, rwidth)

vsini_train = vsini_ol[train_idx]
vsini_cross = vsini_ol[cross_idx]

; Do some learning curves!
data_corrpix = smooth_spec_ol[corrpix,*]

nfeatures_min = 1
nfeatures_max = 15
rdeg = 1

result_sample = {nfeatures:0, $
                 regression_degree:rdeg, $
                 rcoeff:dblarr(nfeatures_max), $
                 rconst:0d0, $
                 vsini_train:vsini_train, $
                 fit_train:vsini_train, $
                 rms_train:0d0, $
                 vsini_cross:vsini_cross, $
                 fit_cross:vsini_cross, $
                 rms_cross:0d0}
                 

niter = nfeatures_max-nfeatures_min + 1

result = replicate(result_sample, niter)
                 

for i = nfeatures_min, nfeatures_max do begin

    ; take just the features we want
    data_mtx = data_corrpix[lindgen(i),*]

    ; Train the data on the training set
    training_set = data_mtx[*,train_idx]
    crossval_set = data_mtx[*,cross_idx]
    
    ; Regress
    rcoeff = regress(training_set, vsini_train, const=rconst, /double)

    ; Apply the fit to all sets
    vsini_fit_train = (training_set ## rcoeff) + rconst
    vsini_fit_cross = (crossval_set ## rcoeff) + rconst
    
    ; Calculate RMS error
    res_train = vsini_train - vsini_fit_train
    res_cross = vsini_cross - vsini_fit_cross
    
    rms_train = sqrt(mean(res_train^2))
    rms_cross = sqrt(mean(res_cross^2))
    
    ; Save the error as a function of # features
    result[i-nfeatures_min].nfeatures = i
    result[i-nfeatures_min].fit_train = vsini_fit_train
    result[i-nfeatures_min].fit_cross = vsini_fit_cross
    result[i-nfeatures_min].rms_train = rms_train
    result[i-nfeatures_min].rms_cross = rms_cross
    result[i-nfeatures_min].rconst = rconst
    
    result[i-nfeatures_min].rcoeff = reform(rcoeff)
    
    
endfor

;;; Plot result
yrmax = max(rms_train) > max(rms_cross)
yrmin = min(rms_train) < min(rms_cross)
yrmax = 2
yrmix = 0


plot, result.nfeatures, result.rms_train, xs=2, yr=[yrmin,yrmax], $
  tit="Learning Curve", $
  xtit="Number of Features", $
  ytit="RMS of Residuals"

oplot, result.nfeatures, result.rms_cross, linest=2

al_legend, ['Training ('+strtrim(ntrain,2)+')', $
            'Cross Val ('+strtrim(ncross,2)+')'], $
  linestyle=[0,2], /right


stop

nvbins = 9

rms = dblarr(n_elements(result), nvbins)
vcmean = dblarr(n_elements(result), nvbins)
vfmean = dblarr(n_elements(result), nvbins)

psopen, 'floor_test.eps', /encapsulated, /color, xs=8, ys=10,/inches
!p.multi=[0,3,5]

for i = 0,nfeatures_max-1 do begin

    resi = result[i]
    
    vcross = resi.vsini_cross
    vfit = resi.fit_cross

    ;diff_i = resi.vsini_cross-resi.fit_cross

    vbin_vec = lindgen(nvbins)+2
    
    bin_idx = value_locate(vbin_vec, resi.vsini_cross)

    for bini = 0, max(bin_idx)-1 do begin
        
        idx = where(bin_idx eq bini, n)
        if n ne 0 then begin
        
            vcross_i = vcross[idx]
            vfit_i = vfit[idx]

            vcmean[i,bini] = mean(vcross_i)
            vfmean[i,bini] = mean(vfit_i)

            rms_i = sqrt(mean((vcross_i-vfit_i)^2))

            rms[i,bini] = rms_i
        endif else begin
            rms[i,bini] = !values.d_nan
        endelse

    endfor
   
    plot, vcmean[i,*], vfmean[i,*], ps=8, xr=[2,10], yr=[0,14], /xs, /ys, tit="RMS for N = "+strtrim(i+1,2)+" features"
    oploterror, vcmean[i,*], vfmean[i,*], rms[i,*], ps=8
    oplot, lindgen(100), lindgen(100), linest=2
    
endfor

psclose
        

;;; Plot data

; nfeat = 6

; idx = nfeat-1

; resstr = result[idx]

; rconst = resstr.rconst

; rcoeff = resstr.rcoeff
; rcoeff = rcoeff[0:idx]

; print, rcoeff
; stop

; rcoeff = transpose(rcoeff)

; full_mtx = data_corrpix[lindgen(nfeat),*]

; full_fit = (full_mtx ## rcoeff) + rconst

; plot, vsini_ol, full_fit, ps=8, xs=2
; oplot, lindgen(100), lindgen(100), linest=2




stop


;;; Other data
;as = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allStar-l31c.2.fits',1)
; alldata_file = '/home/stgilhool/APGBS/spectra/apgbs_datacube.fits'

; ad_full = mrdfits(alldata_file, 1)

; sel_idx = where(ad_full.teff ge t0 and ad_full.teff lt t1, nsel)

; ad = ad_full[sel_idx]

; spectra = ad.spec

; ; Smooth all

;smooth_spec = ml_savgol_spectra(logwl_grid, spectra, width=5)


;apg_sdata = as[ad.allstar_idx]




;;; ID correlated regions

;;; Collect data from correlated regions

;;; Define training set

;;; Run MLR

;;; Check error as a function of CKS vsini (CKS vsini upper limit of
;;; 2km/s for stars with vsini < 1 km/s.  Otherwise, error is ~1km/s)







end
