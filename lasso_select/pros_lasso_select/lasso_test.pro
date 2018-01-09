;+
; NAME: LASSO_TEST
;
;
;
; PURPOSE: Run procedure to test if we can use LASSO regression to
; choose features
;
;
;
; CATEGORY: Analysis
;
;
;
; CALLING SEQUENCE: lasso_test
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

;;;;;;;;;;;;;;;;;;;;
; Function that models Vsini, given design matrix and measured y-values

function lasso_model, params, design_matrix=design_matrix

model = design_matrix ## params

model = reform(model)

return, model

end
;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;
; Function called by mpfit to optimize the model parameters

function lasso_regress, params

common lasso, design_matrix, yvals, err, vis, lasso_opt, lambda
  
; Make the model, given the parameters, x and y data
model = lasso_model(params, design_matrix=design_matrix)

data = yvals

; Calculate residuals
res = data - model

; Calculate deviation
dev = res/err

chi2 = sqrt(total(dev^2, /double))

rmse = sqrt(mean(res^2))

if lasso_opt then begin

    ; Add the penalty
    abs_params = total(abs(params[1:*])) ; not including the constant term
    lasso_penalty = lambda * abs_params

    chi2 = chi2 + lasso_penalty

    ; Make string for plot
    lasso_message = "LASSO is on. Lambda = "+strtrim(lambda,2)+" | abs_params = "+strtrim(abs_params,2)

endif else lasso_message = "LASSO is off"
    
if vis eq 1 then begin

    ; Sort the data by Vsini for visualization
    sorti = sort(data)
    sortd = data[sorti]
    sortm = model[sorti]
    sortr = res[sorti]

    plot, sortd, ps = 8, xtitle = "File number", ytitle = "Vsini", title = lasso_message, /xs
    oplot, sortm, ps = 8, color = !red
    plot, sortr, ps = 8, title = "Residuals with RMSE: "+strtrim(rmse,2), /xs
    wait, 0.001

endif

return, chi2

end

;;;;;;;;;;;;;;;;;;;;


pro lasso_test, LASSO=lasso

common lasso, design_matrix, yvals, err, vis, lasso_opt, lambda

lasso_opt = keyword_set(lasso)

;if n_params() ne 1 then message, "Must input regularization parameter lambda!" 

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

t0 = 4700
t1 = 6600

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
;if n_elements(ad_full[ri])-n_elements(boot_data) eq 48 then print, "N_ELEMENTS reduced by 48, as expected" else print, "Uh oh, overlap overlaps suspect!"
elim_test = cgsetintersection(boot_data.allstar_idx, final_astar_idx, success=elim_check)
if elim_check then print, "Uh oh, there are still overlaps!" else print, "Okay: no stars in the overlap sample present in the boot sample :)"

; Yay!


; Partition
param_arr = [vsini_ol]

nsets = 2L

;divide points between training and cross-val
npts_arr = replicate(nsel_ol/nsets, nsets)
nmod = nsel_ol mod nsets
npts_arr[0] = npts_arr[0] + nmod

;set_str = ml_partition_data(param_arr, npts_arr, /span)

;train_idx = set_str.(0)
;ntrain = n_elements(train_idx) 
;cross_idx = set_str.(1)
;ncross = n_elements(cross_idx) 

;train_idx = lindgen(23)
;cross_idx = lindgen(22) + 23

; Smooth that shit
smooth_spec_ol = ml_savgol_spectra(logwl_grid, spectra_ol, width=5)


; Y-vectors for regression
;vsini_train = vsini_ol[train_idx]
;vsini_cross = vsini_ol[cross_idx]

; Choose features using the old method
; Correlations
corr_vec = ml_get_correlation(smooth_spec_ol, vsini_ol)

; choose regions
rwidth = 9
corrpix = ml_identify_regions(corr_vec, rwidth)

data_corrpix = smooth_spec_ol[corrpix,*]

data = data_corrpix[0:49,*]



nfeatures = n_elements(data[*,0])


 
rdeg = 1 ; just linear regression

; result_sample = {nfeatures:nfeatures, $
;                  lambda:lambda, $
;                  regression_degree:rdeg, $
;                  rcoeff:dblarr(nfeatures), $
;                  rconst:0d0, $
;                  vsini_train:vsini_train, $
;                  fit_train:vsini_train, $
;                  rms_train:0d0, $
;                  vsini_cross:vsini_cross, $
;                  fit_cross:vsini_cross, $
;                  rms_cross:0d0}
                 

;;; Do Multiple Linear Regression

; Add column for constant term
nspec = n_elements(odata) 

if nfeatures gt nspec then message, "Still too many features"

design_matrix = [reform(replicate(1d0,nspec),1,nspec),data]

; Get rid of NaNs
nan_idx = where(finite(design_matrix) eq 0, nancnt)
if nancnt gt 0 then begin
    print, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print, "WARNING: NaNs exist in data matrix. Replacing with zeroes and stopping execution for user review"
    print, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    design_matrix[nan_idx] = 0d0
    stop
endif
    

; get rid of unnecessary stuff
overlaps_data = 0
ad_full = 0
smooth_spec_ol = 0
data = 0

; Add row at bottom to represent the constraint that 
; |p_0| + |p_1| + ... + |p_N-1| le lambda
; Damn, this won't work, actually... I could maybe get around the
; absolute value thing, but I don't know about the inequality... I
; could just do equality, and require that the sum of params ===
; lambda...
; That should actually be okay, OKAY! Now, what about the abs? I
; suppose I could have the design matrix change from iteration to
; iteration. Or I could bias the data and require params to be
; positive? Not sure if/how that would work... Let's try the first way
refit_iteration = 0
refit:

; so let's just add lambda as the last entry in the yvals vector
yvals = vsini_ol

; error
err = replicate(1d0, nspec) ;1 km/s error on cks vsini

; watch output
vis = 1

rcoeff = regress(design_matrix[1:*,*], yvals, const=rconst, /double)

rcoeff = reform(rcoeff)

; Define p0 (guess) and scale for amoeba
nparams = n_elements(design_matrix[*,0])     

guess = [rconst,rcoeff]

scale = 0.5 * guess

ftol = 1d-1
nmax = 200

; If we're on the 2nd (lasso) iteration, use the ordinary regression
; coefficients as guesses
if refit_iteration then begin 

    guess = r
    scale = 3d0*r

    ftol = 1d-3
    nmax = 2000
    
    window, 2

endif else window, 0

!p.multi = [0,1,2]

r = amoeba3(ftol, function_name='lasso_regress', function_value=fval, $
            ncalls=ncalls, nmax=nmax, p0=guess, scale=scale)



!p.multi = 0

if n_elements(r) eq n_elements(guess) then begin

    ; Get the model
    model = lasso_model(r, design_matrix=design_matrix)
    ; in order to get the residuals
    residuals = yvals - model
    ; in order to get RMSE
    rmse = sqrt(mean(residuals^2))

    ; Take chi2 from result
    chi2 = fval[0]

    abs_params = total(abs(r[1:*]))

    if lasso_opt then begin
        ; Check how lasso went
        
        
        lambda_penalty = abs_params * lambda

        lasso_message = "LASSO: Lambda = "+strtrim(lambda,2)+ $
          " | abs_params = "+strtrim(abs_params,2)+" | pen = "+ $
          strtrim(lambda_penalty,2)+" | chi2 = "+strtrim(chi2,2)
        
    endif else lasso_message = "LASSO not used: chi2 = "+strtrim(chi2,2)
        
        ; Plot the result
    window, refit_iteration*2+1
    plot, r[1:*], ps=8, xs=2, ys=2, xtit="Param Number", ytit="Param value", $
      tit="RMSE: "+strtrim(rmse,2)+" | "+lasso_message
    oplot, r[1:*], ps=8, color=!red

    refit_iteration++
    lasso_opt = 1
    lambda = (chi2/abs_params) * 0.5

    
    if refit_iteration eq 1 then goto, refit

endif else begin
    
    print, "AMOEBA failed"
    print, "Result: ", r

endelse

stop

end
