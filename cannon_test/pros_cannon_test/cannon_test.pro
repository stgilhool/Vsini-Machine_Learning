;+
; NAME: CANNON_test
;
;
;
; PURPOSE: Run CANNON-like thing for our approach
;
;
;
; CATEGORY: Analysis
;
;
;
; CALLING SEQUENCE: cannon_test
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


pro cannon_test

common lasso, design_matrix, yvals, err, vis, lasso_opt, lambda

lasso_opt = keyword_set(lasso)

;if n_params() ne 1 then message, "Must input regularization parameter lambda!" 

;;; First, grab some data
;Need spectra, parameters, and map between spectra and parameters
; need cks vsini and apogee spectral types

; Data for the overlaps
overlaps_data = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',1)
logwl_grid = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',0)
; Relevant overlap vectors
teff_ol_cks = overlaps_data.teff_cks
teff_ol_apg = overlaps_data.teff_apg
vsini_ol_cks = overlaps_data.vsini_cks


; Taking just the CKS labels for this test

sel_idx_ol = where(vsini_ol_cks le 25 and vsini_ol_cks ge 1, $
                   nsel_ol)

odata = overlaps_data[sel_idx_ol]


;;; Other data

; Read in APOGEE data cube
alldata_file = '/home/stgilhool/APGBS/spectra/apgbs_datacube.fits'
ad_full = mrdfits(alldata_file, 1)

; Some overlaps are not in data cube due to low S/N. Take those out of
; ol sample
boot_astar_idx = ad_full.allstar_idx
ol_astar_idx = odata.apg_idx

final_astar_idx = cgsetintersection(boot_astar_idx, ol_astar_idx, indices_a=bastar_idx, indices_b=oastar_idx)

; Further pare down the overlap sample
odata = odata[oastar_idx]
teff_ol = odata.teff_cks
teff_ol_apg = odata.teff_apg
vsini_ol = odata.vsini_cks
spectra_ol = odata.spec
error_ol = odata.err
mask_ol = odata.mask

; Smooth that shit, WITH ERRORS!
smooth_err = []
smooth_spec_ol = savgol_custom(logwl_grid, spectra_ol, error_ol, width=5, $
                               savgol_error=smooth_err)

nspec_total = n_elements(odata) 

;;; Data readin finished
print, "Data readin finished"
print, nspec_total, format='("Number of training spectra = ", I0)'
stop
;;;

; Partition into sets. K-fold or leave one out, or two fold for now?
; Just two-fold for now
param_arr = [[vsini_ol], [teff_ol]]

nsets = 2L

;divide points between training and cross-val
npts_arr = replicate(nspec_total/nsets, nsets)
nmod = nspec_total mod nsets
npts_arr[0] = npts_arr[0] + nmod

set_str = ml_partition_data(param_arr, npts_arr, /span)

train_idx = set_str.(0)
ntrain = n_elements(train_idx) 
cross_idx = set_str.(1)
ncross = n_elements(cross_idx) 

;;;;;;;;;;;;;;;;;;;;
;;; TRAIN THE MORTAR
;;;;;;;;;;;;;;;;;;;;



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
    ;scale = 3d0*r
    scale = r

    ftol = 1d-3
    nmax = 5000
    
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
    window, refit_iteration*2+1, xs=1200
    plot, r[1:*], ps=8, xs=2, ys=2, xtit="Param Number", ytit="Param value", $
      tit="RMSE: "+strtrim(rmse,2)+" | "+lasso_message
    oplot, r[1:*], ps=8, color=!red

    lower_idx = where(r[1:*] lt rcoeff, nlow, comp=higher_idx, ncomp=nhigh)

    if nlow gt 0 then begin
        oploterror, lower_idx, rcoeff[lower_idx], $
          rcoeff[lower_idx]-r[lower_idx+1], ps=8, /lobar, /nohat
    endif

    if nhigh gt 0 then begin
        oploterror, higher_idx, rcoeff[higher_idx], $
          r[higher_idx+1]-rcoeff[higher_idx], ps=8, /hibar, /nohat
    endif

    ; Draw dashed line at 0
    oplot, [0,nparams], replicate(0,2), /thick, linest=2

    refit_iteration++
    lasso_opt = 1
    lambda = (chi2/abs_params)

    
    if refit_iteration eq 1 then goto, refit

endif else begin
    
    print, "AMOEBA failed"
    print, "Result: ", r

endelse

stop

end
