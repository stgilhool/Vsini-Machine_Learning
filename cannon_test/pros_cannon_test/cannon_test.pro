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
; Function called by amoeba to optimize s^2 (and consequenty theta)

function training_step, params

common training, label_matrix, dflux, e_dflux, theta_lambda, vis
  
; enforce s_squared is positive
s_squared = abs(params[0])

; At given scatter, calculate the theta
; Scatter is added in quadrature to the flux error term
error_lambda = sqrt(e_dflux^2d0 + s_squared)

theta_fit = regress(label_matrix, dflux, measure_errors=error_lambda, $
                    const=theta_const, yfit=model_dflux, status=rstatus, /double)
; reform column vecs
theta_fit = reform(theta_fit)
model_dflux = reform(model_dflux)

theta_iter = [theta_const, theta_fit]

if rstatus eq 0 then theta_lambda = theta_iter
    
; Calculate the residuals
residuals = model_dflux - dflux

chi2 = total((residuals/error_lambda)^2d0,/double)

; In this formulation, chi2 can be arbitrarily small. I want to
; minimize the difference between chi2/dof and 1

dof = n_elements(dflux)-n_elements(theta_iter) 

chi2_per_dof = chi2/dof

chi2_diff = abs(chi2_per_dof - 1d0)

; penalize chi2 if bad regression
if rstatus ne 0 then chi2_diff = chi2_diff * 1.5

if vis eq 1 then begin

    ; Sort the data by Vsini for visualization
    vsini = reform(label_matrix[1,*])
    ; Sort the data by Teff for visualization
    teff = reform(label_matrix[0,*])

    plot, vsini, dflux, ps=8, symsize=0.5, xs=2, ys=2, xtit="Vsini", ytit="dF/dWL"
    oplot, vsini, model_dflux, ps=8, symsize=0.5, co=!red
    
    plot, teff, dflux, ps=8, symsize=0.5, xs=2, ys=2, xtit="Teff", ytit="dF/dWL"
    oplot, teff, model_dflux, ps=8, symsize=0.5, co=!red

    plot, residuals, ps=8, ytit="Residuals", tit="Scatter = "+strtrim(s_squared, 2)+$
      " | Chi2 = "+strtrim(chi2_per_dof,2), xs=2, ys=2, charsize=2
    oplot, n_elements(residuals)*[-2,2], [0d0,0d0], linest=2, /thick

    wait, 0.001

endif

return, chi2_diff

end

;;;;;;;;;;;;;;;;;;;;


pro cannon_test, VISUALIZE=visualize

common training, label_matrix, dflux, e_dflux, theta_lambda, vis

vis = keyword_set(visualize)
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
error_temp = (error_ol < 10d0)
smooth_spec_ol = savgol_custom(logwl_grid, spectra_ol, error_temp, width=5, $
                               savgol_error=smooth_err)

nspec_total = n_elements(odata) 
npix = n_elements(spectra_ol[*,0]) 
; get rid of unnecessary stuff
overlaps_data = 0
ad_full = 0

;;; Data readin finished
print, "Data readin finished"
print, nspec_total, format='("Number of training spectra = ", I0)'
;stop
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

; Step through each pixel

; If some percentage are masked, skip somehow (inter chip regions
; especially)
; Inter chip: all masked - do something simple (0 slope, high error)
; Bad wavelegth: many masked - probably something similar... take median
;                          of good values, set error high
; Normal wavelength: some masked - process normally. individual pixel
;                                  should have high uncertainty

; Make design matrix
teff_train = teff_ol[train_idx]
vsini_train = vsini_ol[train_idx]

teff_label = teff_train - mean(teff_train)
vsini_label = vsini_train - mean(vsini_train)
constant = replicate(1d0, ntrain)

design_matrix = [transpose(constant), transpose(teff_label), transpose(vsini_label)]
; no constant column for input to regress function
label_matrix = design_matrix[1:*,*]

; Data
slope_data = smooth_spec_ol[*,train_idx]
err_data = smooth_err[*,train_idx]
mask_data = mask_ol[*,train_idx]

; Using all pixel masks for now
mask_bits = lindgen(15)
mask_dec = long(total(2L^mask_bits))


;;; Begin loop through pixels
; Initialize
lambda_mask = bytarr(npix)
theta_arr = dblarr(n_elements(design_matrix[*,0]), npix) 
scatter_vec = dblarr(npix)

for pixnum = 0L, npix-1 do begin

    ; Check if it's the interchip region, or spectacularly bad
    col_mask_dec = reform(mask_data[pixnum,*])
    col_mask_bin = ((mask_dec and col_mask_dec) < 1) ;Set all masked pixels to 1

    n_masked = total(col_mask_bin)

    if n_masked ge (0.9*ntrain) then begin
        ; This pixel is either interchip, or just terrible
        ; Set slope to 0, and err in slope to 1d8
        ; Or maybe just handle at regression step?
        
        print, strtrim(pixnum,2)+" is Bad Pixel, moving on to "+strtrim(pixnum+1,2)
        ; Set scatter to max (don't know what max is yet)
        ; Let's instead flag the pixel
        lambda_mask[pixnum] = 1B
        scatter_vec[pixnum] = !values.d_nan

        ; Set pars to 1,0,0 (1*1 + 0*teff + 0*vsini = 1 (continuum))
        theta_arr[*,pixnum] = [1d0,0d0,0d0]

        continue
    endif


    ;;; Set up amoeba
    ; We are optimizing for theta and s(catter). Actually, just s, at
    ; the amoeba level, I think. At each s, theta is determined by
    ; linear algebra (we can use regress, I think). So we need the
    ; amoeba function to return s, and preferably theta, though we
    ; could just call it with the optimized s to get theta

    
    dflux = reform(slope_data[pixnum,*])
    e_dflux = reform(err_data[pixnum,*])
    
    ftol = 1d-5

    guess = [mean(e_dflux^2, /nan)]  
    scale = guess

    if vis then begin
        window, 0, xs=1400, ys=900, $
          title="Training at Pixel #: "+strtrim(pixnum,2)
        !p.multi = [0,1,3]
    endif

    scatter_lambda = amoeba3(ftol, function_name='training_step', p0=guess, $
                             scale=scale, function_value=fval, $
                             nmax=nmax, ncalls=ncalls)
    
    if vis then !p.multi = 0
    
    ; Test to see if result converged
    if scatter_lambda[0] eq -1 then begin
        
        ; Mask the pixel
        lambda_mask[pixnum] = 2B
        scatter_vec[pixnum] = !values.d_nan

        ; Set pars to 1,0,0 (1*1 + 0*teff + 0*vsini = 1 (continuum))
        theta_arr[*,pixnum] = [1d0,0d0,0d0]

        print, pixnum, format='("WARNING: AMOEBA failed at pixel ",I)'
        print, ''
        ; help, scatter_lambda
;         help, theta_lambda
;         help, guess
;         help, scale
;         print, ''
        
        
        
    endif else begin

        ; Mark pixel as good
        lambda_mask[pixnum] = 0B

        ; Store optimized scatter value
        scatter_i = abs(scatter_lambda[0])
        scatter_vec[pixnum] = scatter_i

        ; Get theta_lambda by calling the function once more
        null_var = training_step(scatter_lambda)
        theta_arr[*,pixnum] = theta_lambda

    endelse

    if ~keyword_set(visualize) then vis = 0

    if pixnum eq 300 or pixnum eq 350 then begin
        
        help, theta_arr
        help, scatter_vec
        help, lambda_mask
        
        vis = 1

        ;stop
    endif

endfor

print, max(scatter_vec, maxidx, /nan), format='("Maximum scatter is: ",F)'

; Set bad pixels to twice maximum scatter

bad_idx = where(lambda_mask ne 0, nbad)

if nbad gt 0 then scatter_vec[bad_idx] = max(scatter_vec, /nan) * 2d0

outstr = {lambda_mask:lambda_mask, scatter_vec:scatter_vec, theta_arr:theta_arr}

print, "Finished"

print, ''

print, "Continue to write structure to disk"

stop

mwrfits, outstr, 'training_param_str_0.fits', /create



stop

;;; Now, the test step
; Using the theta and scatter, we optimize the labels for each
; spectrum.

; With this model, I can trivially retrieve the labels for the test set

slope_data_cross = smooth_spec_ol[*,cross_idx]
err_data_cross = smooth_err[*,cross_idx]
mask_data_cross = mask_ol[*,cross_idx]


;;; Begin loop through pixels
; Initialize
lambda_mask = bytarr(npix)
theta_arr = dblarr(n_elements(design_matrix[*,0]), npix) 
scatter_vec = dblarr(npix)



end
