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

; Now add the row requiring that the sum of the abs values of params =
; lambda
sign_flip_vec = ((params ge 0) * 2) - 1 ;Should be 1 if positive, -1 if negative
; need to zero the first element, since I don't want to include the
; constant offset in the lasso
sign_flip_vec[0] = 0

constraint_row = params * sign_flip_vec

new_design_matrix = [[design_matrix],[constraint_row]]

model = new_design_matrix ## params

model = reform(model)

return, model

end
;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;
; Function called by mpfit to optimize the model parameters

function lasso_regress, params, design_matrix=design_matrix, yvals=yvals, err=err, vis=vis

; Make the model, given the parameters, x and y data
model = lasso_model(params, design_matrix=design_matrix)

data = yvals

; Calculate residuals
res = data - model

; Calculate deviation
dev = res/err

chi2 = sqrt(total(dev^2, /double))

sorti = sort(data)
sortd = data[sorti]
sortm = model[sorti]
sortr = res[sorti]

; Check if lasso constraint is satisfied
param_tot = total(abs(params[1:*]))
lasso_tot = model[-1]

if lasso_tot gt param_tot then begin
    lasso_message = "lasso FAIL"
endif else lasso_message = "lasso COOOOOL"

if vis eq 1 then begin
    
    plot, sortd, ps = 8, xtitle = "File number", ytitle = "Vsini", title = lasso_message, /xs
    oplot, sortm, ps = 8, color = !red
    plot, sortr, ps = 8, title = "Residuals with RMSE: "+strtrim(stddev(res),2), /xs
    wait, 0.001

endif

return, dev

end

;;;;;;;;;;;;;;;;;;;;


pro lasso_test, lambda

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

; Do some learning curves!

; use ALL the pixels! bwahahahah!
; whoops, can't do that. can't have more parameters than spectra
data = smooth_spec_ol[450:500,*]



nfeatures = n_elements(data[*,0])


 
rdeg = 1 ; just linear regression

; first guess lambda
if n_params() eq 0 then lambda = 10

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

data_matrix = [reform(replicate(1d0,nspec),1,nspec),data]

; Get rid of NaNs
nan_idx = where(finite(data_matrix) eq 0, nancnt)
if nancnt gt 0 then begin
    print, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print, "WARNING: NaNs exist in data matrix. Replacing with zeroes and stopping execution for user review"
    print, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    data_matrix[nan_idx] = 0d0
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

; so let's just add lambda as the last entry in the yvals vector
yvals_vsini = vsini_ol
yvals = [yvals_vsini, lambda]

; error
err_data = replicate(1d0, nspec) ;1 km/s error on cks vsini
err = [err_data, 1d-5] ;Make constraint row have tiny error

; watch output
vis = 1

functargs = {design_matrix:data_matrix, $
             yvals:yvals, $
             err:err, $
             vis:vis $
            }

nparams = n_elements(data_matrix[*,0])     
parinfo = replicate({value:0.1d0},nparams)
;parinfo[0].value=20d0

!p.multi = [0,1,2]

r = mpfit('lasso_regress', parinfo=parinfo, functargs=functargs, status=status, dof=dof, bestnorm=chi2, errmsg=errmsg, /quiet)


!p.multi = 0


if status gt 0 then begin

    ; Get the model
    model = lasso_model(r, design_matrix=data_matrix)
    ; in order to get the residuals
    residuals = yvals[0:-2] - model[0:-2]
    ; in order to get RMSE
    rmse = sqrt(mean(residuals^2))


    ; Also, check if lasso constraint is satisfied
    param_tot = total(abs(r[1:*]))
        
    if param_tot gt lambda then begin
        lasso_message = "LASSO constraint violated by "+strtrim(param_tot-lambda,2)+" for lambda = "+strtrim(lambda,2)
    endif else lasso_message = "LASSO constraint satisfied for lambda = "+strtrim(lambda,2)
    
    ; Plot the result
    plot, r[0:-2], ps=8, xs=2, xtit="Param Number", ytit="Param value", tit="RMSE: "+strtrim(rmse,2)+" | "+lasso_message

endif else begin
    
    print, "MPFIT failed with status = "+strtrim(status,2)
    print, errmsg

endelse

stop

end
