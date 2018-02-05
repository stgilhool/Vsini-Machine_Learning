;+
; NAME: CANNON_test_quad
;
;
;
; PURPOSE: Run CANNON-like thing for our approach, quadratic in labels
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
;;; Custom Scatter plot function

pro sg_scatter_plot, x, y, error, COLORS=colors, _EXTRA=ex

plot, x, y, ps=8, _extra=ex, /nodata

foreach color_i, colors, co_idx do begin

    oplot, [x[co_idx]], [y[co_idx]], color=color_i, ps=8
    oploterror, [x[co_idx]], [y[co_idx]], [error[co_idx]], ps=8, errcolor=color_i

endforeach

end

;;;;;;;;;;;;;;;;;;;;
; Function called by amoeba to optimize labels for test spectra

function test_step, params

common test, test_matrix, test_sigma, test_data

nlabels_test = n_elements(params) 

; Loop through to make quadratic terms
quad_terms = []

for label_idx = 0, nlabels_test-1 do begin

    label_i = params[label_idx]

    label_others = params[label_idx:*]

    label_i_terms = label_i * label_others

    quad_terms = [quad_terms, label_i_terms]

endfor

; Full label vector

labels = [1d0, params, quad_terms] 
  
; Do the matrix multiplication (Y = M Labels), where M is test_matrix
; (thetas)
model = test_matrix ## labels

model = reform(model)


residuals = model - test_data

chi2 = total((residuals/test_sigma)^2d0,/double)

; Minimize chi2 (but don't overfit)

dof = n_elements(test_data)-n_elements(labels) 

chi2_per_dof = chi2/dof

chi2_diff = abs(chi2_per_dof - 1d0)



; if vis eq 1 then begin

;     ; Sort the data by Vsini for visualization
;     vsini = reform(label_matrix[1,*])
;     ; Sort the data by Teff for visualization
;     teff = reform(label_matrix[0,*])

;     plot, vsini, dflux, ps=8, symsize=0.5, xs=2, ys=2, xtit="Vsini", ytit="dF/dWL"
;     oplot, vsini, model_dflux, ps=8, symsize=0.5, co=!red
    
;     plot, teff, dflux, ps=8, symsize=0.5, xs=2, ys=2, xtit="Teff", ytit="dF/dWL"
;     oplot, teff, model_dflux, ps=8, symsize=0.5, co=!red

;     plot, residuals, ps=8, ytit="Residuals", tit="Scatter = "+strtrim(s_squared, 2)+$
;       " | Chi2 = "+strtrim(chi2_per_dof,2), xs=2, ys=2, charsize=2
;     oplot, n_elements(residuals)*[-2,2], [0d0,0d0], linest=2, /thick

;     wait, 0.001

; endif

return, chi2_diff

end


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


pro cannon_test_quad, VISUALIZE=visualize, SKIP_OPT=skip_opt, DESCRIPTION=description

common training, label_matrix, dflux, e_dflux, theta_lambda, vis

nlabels_set = 2

vis = keyword_set(visualize)
skip_opt = keyword_set(skip_opt)


; Initialize some stuff
outpath = '/home/stgilhool/Vsini_ML/cannon_test/data_cannon_test/'

log_file = outpath+"test_log.txt"

readcol, log_file, test_num_col, description_col, structure_flag, plot_flag, $
  format='I,A,I,I', delimiter="|", $
  comment="#", /silent, count=ntests

if ntests gt 0 then begin

    ; Last row in test_log
    prev_test_number = test_num_col[-1]
    test_description_default = description_col[-1]
    training_done = structure_flag[-1]
    test_done = plot_flag[-1]

    if training_done and test_done then begin

        test_number_int = prev_test_number + 1
        print, test_number_int, format='("Initiating new test number: ",I)'
        skip_opt = 0

    endif else if training_done and ~test_done then begin
        
        test_number_int = prev_test_number
        print, test_number_int, format='("Continuing test number: ",I)'
        skip_opt = 1

    endif else begin
        
        print, "This wasn't supposed to happen!"
        stop

    endelse


endif else begin

    test_number_int = 0
    skip_opt = 0
    test_description_default = ''
    print, test_number_int, format='("Initiating new test number: ",I)'
    

endelse

terminate = 0

while terminate eq 0 do begin

    desc_opt = ''

    read, desc_opt, $
      prompt="Enter Description ["+test_description_default+"]:"
    
    if desc_opt eq "" then description = test_description_default $
    else description = desc_opt
    
    print, "Description -> "+description

    cont_opt = ''
    
    read, cont_opt, prompt="Continue? (y/n)"
    
    case cont_opt of
        
        ('y' or 'Y' or 'Yes' or 'yes' or 'YES'): begin
            print, "Continuing..."
            print, ""
            terminate = 1
            wait, 2
        end
        
        else: begin
            print, "Re-trying..."
            print, ""
        end
        
    endcase
endwhile


test_number = strtrim(test_number_int,2)

outfile = outpath+'training_param_str_'+test_number+'.fits'


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
feh_ol = odata.feh_cks
logg_ol = odata.logg_cks

; Smooth that shit, WITH ERRORS!
smooth_err = []
error_temp = (error_ol < 1d0)
mskidx = where(mask_ol ne 0)
error_temp[mskidx] = 1d0
spectra_ol[mskidx] = (spectra_ol[mskidx] < 1d0)
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


; train_idx = lindgen(nspec_total)
; cross_idx = lindgen(nspec_total)
; ntrain = nspec_total
; ncross = nspec_total
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
feh_train = feh_ol[train_idx]
logg_train = logg_ol[train_idx]

teff_label = teff_train - mean(teff_train)
vsini_label = vsini_train - mean(vsini_train)
feh_label = feh_train - mean(feh_train)
logg_label = logg_train - mean(logg_train)

constant = replicate(1d0, ntrain)

; Full design matrix

full_matrix_transpose = [[constant],[teff_label],$
                         [vsini_label],[feh_label],[logg_label]]

; Initialize some stuff for making the quadratic terms
quad_matrix_transpose = []
label_string = ['Teff', 'Vsini', 'FeH', 'Logg']

; Loop through to make quadratic terms
for label_idx = 1, nlabels_set do begin

    label_i = full_matrix_transpose[*,label_idx]

    label_others = full_matrix_transpose[*,label_idx:nlabels_set]

    quad_term = rebin(label_i, size(label_others, /dim)) * label_others

    quad_matrix_transpose = [[quad_matrix_transpose],[quad_term]]


    ncrossterms = n_elements(label_others[0,*]) 
    label_i_name = transpose(replicate(label_string[label_idx-1], ncrossterms))
    label_others_name = transpose(label_string[label_idx-1:nlabels_set-1])

    ret_char = transpose(replicate(string(10b),ncrossterms))

    print, [label_i_name, label_others_name, ret_char], $
      format='(5(A," * ",A,:,A))'

endfor


linear_matrix_transpose = full_matrix_transpose[*,0:nlabels_set]

design_matrix_transpose = [[linear_matrix_transpose], [quad_matrix_transpose]]

design_matrix = transpose(design_matrix_transpose)


; no constant column for input to regress function
label_matrix = design_matrix[1:*,*]

nparam = n_elements(design_matrix[*,0]) 
nlabel_columns = n_elements(label_matrix[*,0]) 


; Data
slope_data = smooth_spec_ol[*,train_idx]
err_data = smooth_err[*,train_idx]
mask_data = mask_ol[*,train_idx]

; Using all pixel masks for now
mask_bits = lindgen(15)
mask_dec = long(total(2L^mask_bits))


if skip_opt then goto, skip_training

;;; Begin loop through pixels
; Initialize
lambda_mask = bytarr(npix)
theta_arr = dblarr(nparam, npix) 
scatter_vec = dblarr(npix)

for pixnum = 0L, npix-1 do begin

    ; Check if it's the interchip region, or spectacularly bad
    col_mask_dec = reform(mask_data[pixnum,*])
    col_mask_bin = ((mask_dec and col_mask_dec) < 1) ;Set all masked pixels to 1

    n_masked = total(col_mask_bin)

    if n_masked ge (0.95*ntrain) then begin
        ; This pixel is either interchip, or just terrible
        ; Set slope to 0, and err in slope to 1d8
        ; Or maybe just handle at regression step?
        
        print, strtrim(pixnum,2)+" is Bad Pixel, moving on to "+strtrim(pixnum+1,2)
        ; Set scatter to max (don't know what max is yet)
        ; Let's instead flag the pixel
        lambda_mask[pixnum] = 1B
        scatter_vec[pixnum] = !values.d_nan

        ; Set pars to 1,0,0 (1*1 + 0*teff + 0*vsini = 1 (continuum))
        ;theta_arr[*,pixnum] = [1d0,0d0,0d0]
        theta_arr[*,pixnum] = replicate(0d0, nparam)

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
        theta_arr[*,pixnum] = replicate(0d0, nparam)

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
        if vis then begin
            vis = 0
            null_var = training_step(scatter_lambda)
            vis = 1
        endif else null_var = training_step(scatter_lambda)
        
        theta_arr[*,pixnum] = theta_lambda

    endelse

    if ~keyword_set(visualize) then vis = 0

    if pixnum eq 300 or pixnum eq 350 then begin
        
        help, theta_arr
        help, scatter_vec
        help, lambda_mask
        
        ;vis = 1

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



print, "outfile = "+outfile
print, "test_number = "+test_number
print, ""

;stop

mwrfits, outstr, outfile, /create

openw, lun, log_file, /get_lun, /append

training_done_string = strjoin([test_number,description,'1','0'], '|')

printf, lun, training_done_string
free_lun, lun

;stop

skip_training:

;;; Now, the test step
; Using the theta and scatter, we optimize the labels for each
; spectrum.

; With this model, I can trivially retrieve the labels for the test set

if skip_opt then begin
    training_str = mrdfits(outfile,1)
    lambda_mask = training_str.lambda_mask
    scatter_vec = training_str.scatter_vec
    theta_arr = training_str.theta_arr
endif



slope_data_cross = smooth_spec_ol[*,cross_idx]
err_data_cross = smooth_err[*,cross_idx]
mask_data_cross = mask_ol[*,cross_idx]


;;; Begin loop through spectra
; Initialize

common test, test_matrix, test_sigma, test_data

nlabels_test = nlabels_set

test_status_vec = bytarr(ncross)
test_label_results = dblarr(nlabels_set, ncross) 
test_label_errors = dblarr(nlabels_set, ncross) 

for snum = 0, ncross-1 do begin
    print, "Calculating labels for test spec number: "+strtrim(snum,2)
    ; Grab the spectrum
    test_data = slope_data_cross[*,snum]

    test_err = err_data_cross[*,snum]

    ; FIXME: Should use lambda_mask instead of pixel mask, but the
    ; large scatters should also handle this, so maybe it's fine
    test_mask = mask_data_cross[*,snum]
    test_mask_bin = ((mask_dec and test_mask) < 1)
    test_mask_idx = where(test_mask_bin ne 0, nmask_test, comp=test_unmask_idx, $
                          ncomp=nunmask_test)
    
    ; Jack up error for masked pixels
    ; No, don't do this
    ; test_err_masked = test_err

;     if nmask_test gt 0 then begin
;         test_err_masked[test_mask_idx] = 1d8
;     endif 
    test_err_masked = test_err

    ; Combine errors
    test_sigma = sqrt(test_err_masked^2 + scatter_vec)


    ; Test_matrix (design matrix) is theta_arr
    test_matrix = theta_arr
;    tic
    ; Make covariance matrix for returning label error estimates
    ; It's [ M^T C^-1 M ]^-1
    ;input_covar = diag_matrix(test_sigma^2d0)
    inverse_covar = diag_matrix(1d0/(test_sigma^2d0))
    ;inverse_covar = invert(input_covar, inv_stat, /double)
    
 ;   toc

    dm_trans = transpose(test_matrix)
    
  ;  toc

    ; Do the thing in the brackets
    invco_dm = inverse_covar ## test_matrix
    thing_in_brackets = dm_trans ## invco_dm
    
   ; toc

    ; Get the output covariance matrix
    covar_matrix = invert(thing_in_brackets, covar_stat, /double)

    ;toc

    ;;; We need to use amoeba to optimize the labels
    ; INPUTS: theta_arr, test_sigma, test_data, guess (a vector of
    ; nlabels_set)
    ; FIXME: Add covariance matrix to common block (test) in order to
    ; get label errors back from the test step


    ftol_test = 1d-5
    test_guess_all = [mean(teff_train), mean(vsini_train), mean(feh_train), mean(logg_train)]
    test_scale_all = [1000d0, 10d0, 1d0, 2d0]

    test_guess = test_guess_all[0:nlabels_set-1]
    test_scale = test_scale_all[0:nlabels_set-1]
    nmax_test = 5000


    test_results = amoeba3(ftol_test, function_name='test_step', $
                             p0=test_guess, scale=test_scale, $
                             function_value=fval_test, $
                             nmax=nmax_test, ncalls=ncalls_test)
    


    

    ; Store results
    if n_elements(test_results) eq 1 then begin

        print, "AMOEBA failed in the test step. Results are bad for spec number "+strtrim(snum,2)
        test_label_results[*,snum] = replicate(!values.d_nan, nlabels_set)
        test_label_errors[*,snum] = replicate(!values.d_nan, nlabels_set)
        test_status_vec[snum] = 1

    endif else begin

        test_label_results[*,snum] = test_results

        output_variances = diag_matrix(covar_matrix)
        label_variances = output_variances[1:nlabels_set]
        
        test_label_errors[*,snum] = sqrt(label_variances)
        test_status_vec[snum] = 0

    endelse
endfor


;Teff
teff_cross = teff_ol[cross_idx]
test_teff = reform(test_label_results[0,*]) + mean(teff_train)
test_e_teff = reform(test_label_errors[0,*])


;Vsini
vsini_cross = vsini_ol[cross_idx]
test_vsini = reform(test_label_results[1,*]) + mean(vsini_train)
test_e_vsini = reform(test_label_errors[1,*])


nd_true_idx = where(vsini_cross le 5, nndtrue)
nd_fit_true_idx = where(test_vsini[nd_true_idx] le 5, nndfittrue, comp=nd_fit_false_idx, ncomp=nndfitfalse)

nsuccess = double(nndfittrue)/double(nndtrue)

contamination = double(nndfitfalse)/(double(nndfittrue)+double(nndfitfalse))


if nlabels_set ge 3 then begin
;Feh
feh_cross = feh_ol[cross_idx]
test_feh = reform(test_label_results[2,*]) + mean(feh_train)
test_e_feh = reform(test_label_errors[2,*])
endif


if nlabels_set ge 4 then begin
;Logg
logg_cross = logg_ol[cross_idx]
test_logg = reform(test_label_results[3,*]) + mean(logg_train)
test_e_logg = reform(test_label_errors[3,*])
endif

plotdir = '/home/stgilhool/Vsini_ML/cannon_test/plots_cannon_test/'
plotfile = plotdir + 'cannon_test_'+test_number+'.eps'
;loadct, 34
psopen, plotfile, /encaps, /color, xs=10, ys=8, /inches

!p.multi = [0,nlabels_set, nlabels_set]

for label_num = 0, nlabels_set-1 do begin
    
    case label_num of
        
        0: mag_vec = vsini_cross
        1: mag_vec = teff_cross
        2: mag_vec = feh_cross
        3: mag_vec = logg_cross
        else: message, "Error with case statement"
    endcase

    ;; Vsini
    sg_scatter_plot, vsini_cross, test_vsini-vsini_cross, test_e_vsini, xtit="Vsini (CKS)", ytit="Vsini (Cannon Derivative)", charsize=1.5, colors=bytscl(mag_vec)
    
    oplot, [-100,100], replicate(0d0,2), linest=2, /thick
    
    ;oplot, [-100,100], [-100,100], linest=2, /thick
    ;oploterror, vsini_cross, test_vsini-vsini_cross, test_e_vsini, ps=8
    ;if nfail_test gt 0 then oploterror, vsini_cross[failed_idx], $
    ;test_vsini[failed_idx], test_e_vsini[failed_idx], ps=8, co=!red, errcolor=!red
    
    if nlabels_set ge 2 then begin
                                ;; Teff
        sg_scatter_plot, teff_cross, test_teff-teff_cross, test_e_teff, xtit="Teff (CKS)", ytit="Teff (Cannon Derivative)", charsize=1.5, colors=bytscl(mag_vec)
        
        oplot, [3000,9000], replicate(0d0,2), linest=2, /thick
    endif
    
    if nlabels_set ge 3 then begin
                                ;; Feh
        sg_scatter_plot, feh_cross, test_feh-feh_cross, test_e_feh, xtit="[Fe/H] (CKS)", ytit="[Fe/H] (Cannon Derivative)", charsize=1.5, colors=bytscl(mag_vec)
        
        oplot, [-3,3], replicate(0d0,2), linest=2, /thick
    endif
    
    if nlabels_set ge 4 then begin
                                ;; Logg
        sg_scatter_plot, logg_cross, test_logg-logg_cross, test_e_logg, xtit="Logg (CKS)", ytit="Logg (Cannon Derivative)", charsize=1.5, colors=bytscl(mag_vec)
        
        oplot, [2,7], replicate(0d0,2), linest=2, /thick
    endif
    
endfor

xyouts, 0.5, 0.95, "Vsini Non-detection Contamination = "+strtrim(sigfig(contamination,3),2),/alignment, /normal

psclose

; Write to the log file that the plot is made
openw, lun, log_file, /get_lun, /append

; get position of pointer at end of file
point_lun, -1*lun, eof_pos

; move pointer back two characters (carriage return and final '0')
point_pos = eof_pos - 2
point_lun, lun, point_pos

; Delete the final 0
truncate_lun, lun

; Write the final 1
printf, lun, '1'
free_lun, lun

;stop


end


; Mark failed guys in red
; failed_idx = where(test_status_vec ne 0, nfail_test)


; ;window, 1, tit="RESULTS!!!", xs=1500, ys=900
; plotdir = '/home/stgilhool/Vsini_ML/cannon_test/plots_cannon_test/'
; plotfile = plotdir + 'cannon_test_'+test_number+'.eps'
; psopen, plotfile, /encaps, /color, xs=10, ys=8, /inches

; !p.multi = [0,1,2]

; plot, vsini_cross, test_vsini, ps=8, $
;   tit="Quadratic in "+strtrim(nlabels_set,2)+" labels", $
;   xtit="Vsini (CKS)", ytit="Vsini (Cannon Derivative)", charsize=1.5
; oplot, [-100,100], [-100,100], linest=2, /thick
; oploterror, vsini_cross, test_vsini, test_e_vsini, ps=8
; if nfail_test gt 0 then oploterror, vsini_cross[failed_idx], $
;   test_vsini[failed_idx], test_e_vsini[failed_idx], ps=8, co=!red, errcolor=!red


; plot, teff_cross, test_teff, ps=8, xtit="Teff (CKS)", ytit="Teff (Cannon Derivative)", charsize=1.5, /yno
; oplot, [-1d5, 1d5], [-1d5,1d5], linest=2, /thick
; oploterror, teff_cross, test_teff, test_e_teff, ps=8
; if nfail_test gt 0 then oploterror, teff_cross[failed_idx], $
;   test_teff[failed_idx], test_e_teff[failed_idx], ps=8, co=!red, errcolor=!red
