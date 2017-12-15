; First attempt to apply ML technique to all ~35000 spectra ('halfset'
; is antiquated)
pro apgbs_halfset_0

overlaps_data = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',1)
logwl_grid = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',0)

;;; Select Full Data Set
; Select stars of a given temperature
cks_teff = overlaps_data.teff_cks
cks_vsini = overlaps_data.vsini_cks


; Get temperatures that define the bin
teff_bin_0 = 4600L
teff_bin_1 = 7000L

; Select the spectra
sel_idx = where(cks_teff ge teff_bin_0 and $
                cks_teff lt teff_bin_1 and $
                overlaps_data.overlap_idx ne 323 and $
                cks_vsini ge 4 and $
                cks_vsini le 20, nsel)

; Truncate the data structures
data = overlaps_data[sel_idx]
spectra = data.spec

; Get the parameter vectors
vsini = data.vsini_cks
teff = data.teff_cks
feh = data.feh_cks
;;;;;;;;;;;;;;;;;;;;


;;; Get derivatives, compute correlations, and choose regions

; get smoothed derivative
smooth_spectra = ml_savgol_spectra(logwl_grid, spectra, width=5)

; compute correlations
corr_vec = ml_get_correlation(smooth_spectra, vsini)

; choose regions
width = 5
corrpix = ml_identify_regions(corr_vec, width)
;;;;;;;;;;;;;;;;;;;;



; training set
training_set_idx = lindgen(nsel)

vsini_train = vsini
    
data_corrpix = smooth_spectra[corrpix,*]

nfeatures = 15                 
feat_idx = lindgen(nfeatures)
;feat_idx = [2,3]

; take just the features we want
data_mtx = data_corrpix[feat_idx,*]

; Train the data on the training set
;training_set = data_mtx[*,training_set_idx]
training_set = data_mtx
    ;crossval_set = data_mtx[*,crossval_set_idx]
    ;test_set = data_mtx[*,test_set_idx]

    ; Do 2nd degree
    training_set = [training_set, training_set^2d0]
    ;crossval_set = [crossval_set, crossval_set^2d0]
    ;test_set = [test_set, test_set^2d0]

rcoeff = regress(training_set, vsini_train, const=rconst, /double)

    ; Apply the fit to all sets
vsini_fit_train = (training_set ## rcoeff) + rconst
    ;vsini_fit_cross = (crossval_set ## rcoeff) + rconst
    ;vsini_fit_test = (test_set ## rcoeff) + rconst

vsini_fit_train = reform(vsini_fit_train)
    ;vsini_fit_cross = reform(vsini_fit_cross)
    ;vsini_fit_test = reform(vsini_fit_test)


    ; Calculate RMS error
res_train = vsini_train - vsini_fit_train
    ;res_cross = vsini_cross - vsini_fit_cross
    ;res_test = vsini_test - vsini_fit_test

rms_train = sqrt(mean(res_train^2))
    ;rms_cross = sqrt(mean(res_cross^2))
    ;rms_test = sqrt(mean(res_test^2))


;;;; Read in big thing
alldata_file = '/home/stgilhool/APGBS/spectra/apgbs_datacube.fits'
ad = mrdfits(alldata_file, 1)

; apply to alldata_file
smooth_data = ml_savgol_spectra(logwl_grid, ad.spec, width=5)

data_corrpix_all = smooth_data[corrpix[feat_idx], *]

data_set = [data_corrpix_all, data_corrpix_all^2d0]
;data_set = data_corrpix_all

vsini_result = (data_set ## rcoeff) + rconst

plot, ad.vsini_apg, vsini_result, ps=8
oplot, lindgen(1000), lindgen(1000), /thick, linest=2
stop
window, 0, xs=1200, ys=800
plot, ad.vsini_apg, vsini_result, ps=8, yr=[0,150], xr=[0,150], $
  xtit = "Vsini from ASPCAP pipeline", ytit="Vsini from CKS fit", charsize=2
oplot, lindgen(1000), lindgen(1000), /thick, linest=2
stop
plot, ad.vsini_apg, vsini_result, ps=8, yr=[5,20], xr=[5,20], $
  xtit = "Vsini from ASPCAP pipeline", ytit="Vsini from CKS fit", charsize=2
oplot, lindgen(1000), lindgen(1000), /thick, linest=2
stop
plot, vsini_train, vsini_fit_train, ps=8, xtit="Vsini (training)", ytit="Vsini_fit (for training data)", xs=2, ys=2
oplot, lindgen(1000), lindgen(1000), /thick, linest=2
stop

feati = 0
feat_view:

plot, vsini_train, training_set[feati,*], ps=8, xtit="Vsini (training)", ytit="X_"+strtrim(feati,2)
feati = ''
read, feati, prompt="% Which chunk? [0-"+strtrim(n_elements(training_set[*,0]) ,2)+"]"

if ((feati ge 0) && (feati lt n_elements(training_set[*,0]))) then goto, feat_view

stop
        

end
