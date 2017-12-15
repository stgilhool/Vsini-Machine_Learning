pro plot_learning_curves


cks_file_data = mrdfits('/home/stgilhool/CKS/spectra/allspec.fits', 1)
cks_stellar_data = mrdfits('/home/stgilhool/CKS/spectra/allspec.fits', 2)
spec_arr = mrdfits('/home/stgilhool/CKS/spectra/allspec.fits', 3)
wl_grid = mrdfits('/home/stgilhool/CKS/spectra/allspec.fits', 4)

;;; Select Full Data Set
; Select stars of a given temperature
cks_teff = cks_stellar_data.teff
cks_vsini = cks_stellar_data.vsini
cks_logg = cks_stellar_data.logg

; Get temperatures that define the bin
teff_bin_0 = 4600L
teff_bin_1 = 5000L

; Select the spectra
sel_idx = where(cks_teff ge teff_bin_0 and $
                cks_teff lt teff_bin_1 and $
                cks_logg ge 3.5 and $
                cks_vsini le 85, nsel)

; Truncate the data structures
file_info = cks_file_data[sel_idx]
star_info = cks_stellar_data[sel_idx]
spectra = spec_arr[*,sel_idx]

; Get the parameter vectors
vsini = star_info.vsini
teff = star_info.teff
logg = star_info.logg
feh = star_info.feh
;;;;;;;;;;;;;;;;;;;;


;;; Partition Data
param_arr = [[vsini], $
             [teff], $
             [logg], $
             [feh]]

nsets = 3L

; divide all points evenly among the sets
npts_arr = replicate(nsel/nsets, nsets)
; add any remaining points to the first set
nmod = nsel mod nsets
npts_arr[0] = npts_arr[0] + nmod

; OVERWRITING
ntrain = 66L
ncross = 20L
ntest = 20L
npts_arr = [ntrain, ncross, ntest]

set_str = ml_partition_data(param_arr, npts_arr, /span)
;;;;;;;;;;;;;;;;;;;;


;;; Get derivatives, compute correlations, and choose regions

; get smoothed derivative
smooth_spectra = ml_savgol_spectra(wl_grid, spectra)

; compute correlations
corr_vec = ml_get_correlation(smooth_spectra, vsini)

; choose regions
width = 9
corrpix = ml_identify_regions(corr_vec, width)
;;;;;;;;;;;;;;;;;;;;



;;; Now, do regression, get rms for each set as a function of #
;;; features


    ; training set
    training_set_idx = set_str.(0)

    ; cross-val set
    crossval_set_idx = set_str.(1)

    ; test set
    test_set_idx = set_str.(2)


    vsini_train = vsini[training_set_idx]
    vsini_cross = vsini[crossval_set_idx]
    vsini_test = vsini[test_set_idx]



data_corrpix = smooth_spectra[corrpix,*]

nfeatures_min = 1
nfeatures_max = 15

result_sample = {nfeatures:0, $
                 regression_degree:1, $
                 vsini_train:vsini_train, $
                 fit_train:vsini_train, $
                 rms_train:0d0, $
                 vsini_cross:vsini_cross, $
                 fit_cross:vsini_cross, $
                 rms_cross:0d0, $
                 vsini_test:vsini_test, $
                 fit_test:vsini_test, $
                 rms_test:0d0}

niter = nfeatures_max-nfeatures_min + 1

result = replicate(result_sample, niter)
                 

for i = nfeatures_min, nfeatures_max do begin

    ; take just the features we want
    data_mtx = data_corrpix[lindgen(i),*]

    ; Train the data on the training set
    training_set = data_mtx[*,training_set_idx]
    crossval_set = data_mtx[*,crossval_set_idx]
    test_set = data_mtx[*,test_set_idx]

    ; Do 2nd degree
    training_set = [training_set, training_set^2d0]
    crossval_set = [crossval_set, crossval_set^2d0]
    test_set = [test_set, test_set^2d0]

    rcoeff = regress(training_set, vsini_train, const=rconst, /double)

    ; Apply the fit to all sets
    vsini_fit_train = (training_set ## rcoeff) + rconst
    vsini_fit_cross = (crossval_set ## rcoeff) + rconst
    vsini_fit_test = (test_set ## rcoeff) + rconst

    ; Calculate RMS error
    res_train = vsini_train - vsini_fit_train
    res_cross = vsini_cross - vsini_fit_cross
    res_test = vsini_test - vsini_fit_test

    rms_train = sqrt(mean(res_train^2))
    rms_cross = sqrt(mean(res_cross^2))
    rms_test = sqrt(mean(res_test^2))

    ; Save the error as a function of # features
    result[i-nfeatures_min].nfeatures = i
    result[i-nfeatures_min].fit_train = vsini_fit_train
    result[i-nfeatures_min].fit_cross = vsini_fit_cross
    result[i-nfeatures_min].fit_test = vsini_fit_test
    result[i-nfeatures_min].rms_train = rms_train
    result[i-nfeatures_min].rms_cross = rms_cross
    result[i-nfeatures_min].rms_test = rms_test
    
endfor


;;; Plot result
yrmax = max(rms_train) > max(rms_cross) > max(rms_test)
yrmin = min(rms_train) < min(rms_cross) < min(rms_test)

psopen, '/home/stgilhool/CKS/plots/learning_curves_moretraining_semiquad.eps', /encaps, $
  /color, xs=8, ys=8, /inches

plot, result.nfeatures, result.rms_train, xs=2, yr=[yrmin,yrmax], $
  tit="Learning Curve for 'Quadratic' Fit", $
  xtit="Number of Features", $
  ytit="RMS of Residuals"

oplot, result.nfeatures, result.rms_cross, linest=2
oplot, result.nfeatures, result.rms_test, linest=3

al_legend, ['Training ('+strtrim(ntrain,2)+')', $
            'Cross Val ('+strtrim(ncross,2)+')', $
            'Test ('+strtrim(ntest,2)+')'], $
  linestyle=[0,2,3], /right

psclose


stop

end
