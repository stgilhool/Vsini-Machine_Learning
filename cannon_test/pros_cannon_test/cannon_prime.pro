;+
; NAME: Cannon_prime
;
;
;
; PURPOSE: Measure Vsini for APOGEE spectrum with CANNON-like approach
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
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
; Function called by amoeba to optimize labels for test spectra

function test_step, params

common test, test_matrix, test_sigma, test_data, model_type

nlabels_test = n_elements(params) 

; Loop through to make quadratic terms, if necessary
if model_type eq 'quadratic' then begin
    quad_terms = []
    
    for label_idx = 0, nlabels_test-1 do begin
        
        label_i = params[label_idx]
        
        label_others = params[label_idx:*]
        
        label_i_terms = label_i * label_others
        
        quad_terms = [quad_terms, label_i_terms]
        
    endfor

    ; Full label vector

    labels = [1d0, params, quad_terms] 
    
endif else labels = [1d0, params]
  
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


return, chi2_diff

end

;;;;;;;;;;;;;;;;;;;;
; Function called by amoeba to optimize s^2 (and consequenty theta)

function training_step, params

common training, label_matrix, flux, e_flux, theta_lambda, vis
  
; enforce s_squared is positive
s_squared = abs(params[0])

; At given scatter, calculate the theta
; Scatter is added in quadrature to the flux error term
error_lambda = sqrt(e_flux^2d0 + s_squared)

theta_fit = regress(label_matrix, flux, measure_errors=error_lambda, $
                    const=theta_const, yfit=model_flux, status=rstatus, /double)
; reform column vecs
theta_fit = reform(theta_fit)
model_flux = reform(model_flux)

theta_iter = [theta_const, theta_fit]

if rstatus eq 0 then theta_lambda = theta_iter
    
; Calculate the residuals
residuals = model_flux - flux

chi2 = total((residuals/error_lambda)^2d0,/double)

; In this formulation, chi2 can be arbitrarily small. I want to
; minimize the difference between chi2/dof and 1

dof = n_elements(flux)-n_elements(theta_iter) 

chi2_per_dof = chi2/dof

chi2_diff = abs(chi2_per_dof - 1d0)

; penalize chi2 if bad regression
if rstatus ne 0 then chi2_diff = chi2_diff * 1.5

if vis eq 1 then begin

    ; Sort the data by Vsini for visualization
    ;vsini = reform(label_matrix[1,*])
    ; Sort the data by Teff for visualization
    ;teff = reform(label_matrix[0,*])

    ;plot, vsini, dflux, ps=8, symsize=0.5, xs=2, ys=2, xtit="Vsini", ytit="dF/dWL"
    ;oplot, vsini, model_dflux, ps=8, symsize=0.5, co=!red
    
    ;plot, teff, dflux, ps=8, symsize=0.5, xs=2, ys=2, xtit="Teff", ytit="dF/dWL"
    ;oplot, teff, model_dflux, ps=8, symsize=0.5, co=!red

    plot, residuals, ps=8, ytit="Residuals", tit="Scatter = "+ $
      strtrim(s_squared, 2)+$
      " | Chi2 = "+strtrim(chi2_per_dof,2), xs=2, ys=2, charsize=2
    oplot, n_elements(residuals)*[-2,2], [0d0,0d0], linest=2, /thick

    wait, 0.001

endif

return, chi2_diff

end

;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;
;;; MAIN
;;;;;;;;;;;;;;;;;;;;


pro cannon_prime, VISUALIZE=visualize

;;;INITIALIZE
;; Get Settings from input structure or input file
;; Set settings interactively, if no input structure or file
;; Allow just the writing of settings?
;; Open and write a log message
;;
;; Settings:
;; MODEL
;; ; which labels/number of labels
;; ; model type (linear function of some combination of labels)
;; DATA
;; ; apStar or aspcapStar
;; ; normalization method
;; ; which spectra/which cut criteria
;; ; spectra
;; ; flux errors
;; ; mask
;; ; 4 cks labels



common training, label_matrix, flux, e_flux, theta_lambda, vis

; Hardcode the settings
label_names    = ['Vsini', 'Teff', '[Fe/H]', 'logg']
shortnames     = [   'v',     't',      'f',    'g']
label_settings = [     1,        1,      0,       0]

model_names    = ['linear', 'quadratic']
model_settings = [       1,           0]

if total(model_settings) ne 1 then message, "One and only one model must be chosen"


nlabels_set = total(label_settings)

vis = keyword_set(visualize)

; Initialize some stuff, and ask user for description
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
        
        ('n' or 'N' or 'No' or 'no' or 'NO'): begin
            print, "Terminating..."
            print, ""
            return
        end


        else: begin
            print, "Re-trying..."
            print, ""
        end
        
    endcase
endwhile


test_number = strtrim(test_number_int,2)

outfile = outpath+'training_param_str_'+test_number+'.fits'
median_outfile = outpath+'training_median_str_'+test_number+'.fits'
;;; First, grab some data
;Need spectra, parameters, and map between spectra and parameters

; Data for the overlaps
overlaps_data = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',1)
logwl_grid = mrdfits('/home/stgilhool/APGBS/spectra/overlaps/apgbs_overlaps_data.fits',0)
; Relevant overlap vectors
vsini_ol_cks = overlaps_data.vsini_cks


; Taking just the CKS labels for this test

;sel_idx_ol = where(vsini_ol_cks le 2, $
;                   nsel_ol)
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

; How about we try the CANNON without any fake spectra?

; Smooth that shit, WITH ERRORS!
smooth_err = []
max_flux_error = 1d0
error_temp = (error_ol < max_flux_error)
mskidx = where(mask_ol ne 0)
error_temp[mskidx] = max_flux_error
spectra_ol[mskidx] = (spectra_ol[mskidx] < 1d0)
smooth_spec_ol = savgol_custom(logwl_grid, spectra_ol, error_temp, width=5, degree=4, $
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


;;; LEAVE-ONE-OUT
ncross = 1L
ntrain = nspec_total - 1L

spec_idx_full = lindgen(nspec_total)

;;;;;;;;;;;;;;;;;;;;
;;; TRAIN THE CANNON
;;;;;;;;;;;;;;;;;;;;

;;; Do Multiple Linear Regression

; Step through each pixel

; Inter chip: all masked - do something simple (0 slope, high error)
; Bad wavelegth: many masked - probably something similar... take median
;                          of good values, set error high
; Normal wavelength: some masked - process normally. individual pixel
;                                  should have high uncertainty

; FIXME : stuff copied here that I can't skip in the debug step
    ; Using all pixel masks for now
mask_bits = lindgen(15)
mask_dec = long(total(2L^mask_bits))
label_set_idx = where(label_settings eq 1)

;goto, debug

;;; Center and standardize label data
vsini_mean = mean(vsini_ol)
vsini_sig = stddev(vsini_ol)
vsini_all = (vsini_ol-vsini_mean)/(2d0*vsini_sig)

teff_mean = mean(teff_ol)
teff_sig = stddev(teff_ol)
teff_all = (teff_ol-teff_mean)/(2d0*teff_sig)

feh_mean = mean(feh_ol)
feh_sig = stddev(feh_ol)
feh_all = (feh_ol-feh_mean)/(2d0*feh_sig)

logg_mean = mean(logg_ol)
logg_sig = stddev(logg_ol)
logg_all = (logg_ol-logg_mean)/(2d0*logg_sig)

max_loop_iter = nspec_total-1

 if skip_opt then begin
    
     ; Make full design matrix if skipping training setp
     constant = replicate(1d0, nspec_total)
    
     ; Contruct full design matrix
    
     full_matrix_transpose = [[constant],[vsini_all], [teff_all],$
                              [feh_all],[logg_all]]
    
     label_set_idx = where(label_settings eq 1)
    
     linear_matrix_transpose = full_matrix_transpose[*,[0,label_set_idx+1]]
    
    
     ; Add quadratic terms, if necessary
     if model_settings[1] eq 1 then begin
         ; Initialize some stuff for making the quadratic terms
         quad_matrix_transpose = []
        
         ; Loop through to make quadratic terms
         for label_idx = 1, nlabels_set do begin
            
             label_i = full_matrix_transpose[*,label_idx]
            
             label_others = full_matrix_transpose[*,label_idx:nlabels_set]
            
             quad_term = rebin(label_i, size(label_others, /dim)) * label_others
            
             quad_matrix_transpose = [[quad_matrix_transpose],[quad_term]]
            
         endfor
        
         design_matrix_transpose = [[linear_matrix_transpose], [quad_matrix_transpose]]
        
     endif else if model_settings[0] eq 1 then begin
        
         ; linear case
         design_matrix_transpose = linear_matrix_transpose
        
     endif else message, "Model incorrectly specified"
    
     design_matrix = transpose(design_matrix_transpose)
    
     ; no constant column for input to regress function
     label_matrix = design_matrix[1:*,*]
    
     nparam = n_elements(design_matrix[*,0]) 
     nlabel_columns = n_elements(label_matrix[*,0]) 
    
    
     ; Data
 ;    flux_data = spectra_ol[*,train_idx]
     flux_data = smooth_spec_ol
     ;err_data = error_ol[*,train_idx]
     err_data = smooth_err
     mask_data = mask_ol
    
     ; Using all pixel masks for now
     mask_bits = lindgen(15)
     mask_dec = long(total(2L^mask_bits))
    
     goto, skip_training

 endif    
    
for iter = 0, nspec_total-1 do begin

;for iter = 61, max_loop_iter do begin

    ; Shift the index matrix by iter entries
    spec_idx_shift = shift(spec_idx_full, iter)
    
    ; 0th entry is the cross-val index, the rest are training
    cross_idx = spec_idx_shift[0]
    train_idx = spec_idx_shift[1:*]

    
    ; Make design matrix
    teff_train = teff_all[train_idx]
    vsini_train = vsini_all[train_idx]
    feh_train = feh_all[train_idx]
    logg_train = logg_all[train_idx]
    
    teff_label = teff_train 
    vsini_label = vsini_train 
    feh_label = feh_train 
    logg_label = logg_train 
    
    constant = replicate(1d0, ntrain)
    
    ; Contruct full design matrix
    
    full_matrix_transpose = [[constant],[vsini_label], [teff_label],$
                             [feh_label],[logg_label]]
    
    label_set_idx = where(label_settings eq 1)
    
    linear_matrix_transpose = full_matrix_transpose[*,[0,label_set_idx+1]]
    
    
    ; Add quadratic terms, if necessary
    if model_settings[1] eq 1 then begin
        ; Initialize some stuff for making the quadratic terms
        quad_matrix_transpose = []
        
        ; Loop through to make quadratic terms
        for label_idx = 1, nlabels_set do begin
            
            label_i = full_matrix_transpose[*,label_idx]
            
            label_others = full_matrix_transpose[*,label_idx:nlabels_set]
            
            quad_term = rebin(label_i, size(label_others, /dim)) * label_others
            
            quad_matrix_transpose = [[quad_matrix_transpose],[quad_term]]
            
        endfor
        
        design_matrix_transpose = [[linear_matrix_transpose], [quad_matrix_transpose]]
        
    endif else if model_settings[0] eq 1 then begin
        
        ; linear case
        design_matrix_transpose = linear_matrix_transpose
        
    endif else message, "Model incorrectly specified"
    
    design_matrix = transpose(design_matrix_transpose)
    
    ; no constant column for input to regress function
    label_matrix = design_matrix[1:*,*]
    
    nparam = n_elements(design_matrix[*,0]) 
    nlabel_columns = n_elements(label_matrix[*,0]) 
    
    
    ; Data
;    flux_data = spectra_ol[*,train_idx]
    flux_data = smooth_spec_ol[*,train_idx]
    ;err_data = error_ol[*,train_idx]
    err_data = smooth_err[*,train_idx]
    mask_data = mask_ol[*,train_idx]
    
    ; Using all pixel masks for now
    mask_bits = lindgen(15)
    mask_dec = long(total(2L^mask_bits))
    
    
    ;if skip_opt then goto, skip_training
    
    ;;; Begin loop through pixels
    ; Initialize
    lambda_mask = bytarr(npix)
    theta_arr = dblarr(nparam, npix) 
    scatter_vec = dblarr(npix)
    
    for pixnum = 0L, npix-1 do begin
        
        ; Check if it's the interchip region, or spectacularly bad
        col_mask_dec = reform(mask_data[pixnum,*])
        col_mask_bin = ((mask_dec and col_mask_dec) < 1) ;Set all masked pixels to 1
        
        col_mask_idx = where(col_mask_bin eq 1, n_masked, comp=gd_idx, ncomp=ngood)
        
        if n_masked ge (floor(0.9*ntrain)) then begin
            ; This pixel is either interchip, or just terrible
            
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
        
        
        flux = reform(flux_data[pixnum,*])
        e_flux = reform(err_data[pixnum,*])
        
        ; Set flagged pixel error high
        ;e_flux[col_mask_idx] = 1d8
        ftol = 1d-5
        
        guess = [median(e_flux^2)]  
        scale = guess
        
        if vis then begin
            window, 0, xs=1400, ys=900, $
              title="Training at Pixel #: "+strtrim(pixnum,2)
            ;!p.multi = [0,1,3]
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
            
            ;print, pixnum, format='("WARNING: AMOEBA failed at pixel ",I)'
            ;print, ''
            
            
        endif else begin
            
            ; Store optimized scatter value
            scatter_i = abs(scatter_lambda[0])
            
            
            ; Get theta_lambda by calling the function once more
            if vis then begin
                vis = 0
                null_var = training_step(scatter_lambda)
                vis = 1
            endif else null_var = training_step(scatter_lambda)
            
            if scatter_i eq 0 then begin
                
                lambda_mask[pixnum] = 3B
                scatter_vec[pixnum] = !values.d_nan
                theta_arr[*,pixnum] = replicate(0d0, nparam)
                
              ;  print, pixnum, format='("WARNING: Scatter is 0 at pixel ",I)'
               ; print, ''
                
            endif else if total(finite(theta_lambda)) ne nparam then begin

                lambda_mask[pixnum] = 4B
                scatter_vec[pixnum] = !values.d_nan
                theta_arr[*,pixnum] = replicate(0d0, nparam)
                
                ;print, pixnum, format='("WARNING: parameters undefined at pixel ",I)'
                ;print, ''
                ;wait, 0.5
                
            endif else begin

                ; Mark pixel as good
                lambda_mask[pixnum] = 0B
                
                scatter_vec[pixnum] = scatter_i
                theta_arr[*,pixnum] = theta_lambda
                
            endelse
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
    
    print, max(scatter_vec, /nan), format='("Maximum scatter is: ",D0)'
    
    ; Set bad pixels to twice maximum scatter
    
    bad_idx = where(lambda_mask ne 0, nbad)
    
    if nbad gt 0 then scatter_vec[bad_idx] = max(scatter_vec, /nan) * 2d0
    
    outstr = {training_idx:train_idx, $
              cross_idx:cross_idx, $
              lambda_mask:lambda_mask, $
              scatter_vec:scatter_vec, $
              theta_arr:theta_arr}
    
     ;if iter eq 0 then outstr_loo = replicate(outstr, nspec_total) $
     ;else outstr_loo[iter] = outstr
    if iter eq 0 then outstr_loo = replicate(outstr, max_loop_iter+1) $
    else outstr_loo[iter] = outstr
    ; if iter eq 61 then begin
;         outstr_loo_in = mrdfits('/home/stgilhool/Vsini_ML/cannon_test/data_cannon_test/training_param_str_66.fits',1)
;         outstr_loo = replicate(outstr, max_loop_iter+1)
;         outstr_loo[0] = outstr_loo_in
;     endif else outstr_loo[iter] = outstr

    print, "Finished for iteration number: "+ $
      strtrim(iter+1,2)+"/"+strtrim(nspec_total,2)
    
    print, ''
    
    print, "Continue to write structure to disk"
    
    print, "outfile = "+outfile
    print, "test_number = "+test_number
    print, ""
    
    ; Write results, as we go
    if iter eq 0 then begin
        mwrfits, outstr, outfile, /create
    endif else mwrfits, outstr, outfile
    
    ; Write to log file, and overwrite results with full results, if finished
    ;if iter eq nspec_total-1 then begin
    if iter eq max_loop_iter then begin

        ; Overwrite with array of structures
        mwrfits, outstr_loo, outfile, /create
        
        ; Log
        openw, lun, log_file, /get_lun, /append
    
        training_done_string = strjoin([test_number,description,'1','0'], '|')
        
        printf, lun, training_done_string
        free_lun, lun
    endif
    
endfor

;debug:
skip_training:

outstr_loo = mrdfits(outfile, 1)


; Get median for all params and pixels
theta_arr_loo = outstr_loo.theta_arr
lambda_mask_loo = outstr_loo.lambda_mask
scatter_vec_loo = outstr_loo.scatter_vec

; help, theta_arr_loo
; help, lambda_mask_loo
; help, scatter_vec_loo


; Nan bad ones
 for pixnum = 0, npix-1 do begin

     mask_spec = where(reform(lambda_mask_loo[pixnum,*]) ne 0, nmask, comp=gspec)

     if nmask gt max_loop_iter/2 then begin

         ; lambda_mask_loo[pixnum,mask_spec] = 5B
;          scatter_vec_loo[pixnum,mask_spec] = !values.d_nan
;          theta_arr_loo[*,pixnum,mask_spec] = !values.d_nan

          lambda_mask_loo[pixnum,*] = 5B
          scatter_vec_loo[pixnum,*] = !values.d_nan
          theta_arr_loo[*,pixnum,*] = !values.d_nan

     endif

 endfor

theta_arr_median = median(theta_arr_loo, dim=3)
scatter_vec_median = median(scatter_vec_loo, dim=2)

;;; ABORTED ATTEMPT TO DO WEIGHTED AVERAGE
; weights = 1d0/scatter_vec_loo
; weights_tmp = reform(weights, [1, size(weights,/dim)])
; weights_arr = rebin(weights_tmp, size(theta_arr, /dim)]

; theta_arr_weighted = weights_arr * theta_arr_weighted

; weights_tot_arr = total(weights_arr, 3, /double)
; theta_arr_wtot = total(theta_arr_weighted, 3, /double)

; theta_arr_median = theta_arr_wtot/weights_tot_arr

max_scatter = max(scatter_vec_median, /nan)

lambda_mask_idx = where((~finite(scatter_vec_median)) or (~finite(theta_arr_median[0,*])), nmask, comp=good_idx)
lambda_mask_final = replicate(0B, npix)


if nmask ne 0 then begin
    
    theta_arr_median[*,lambda_mask_idx] = 0d0
    scatter_vec_median[lambda_mask_idx] = max_scatter * 2d0
    lambda_mask_final[lambda_mask_idx] = 1B

endif

theta_arr_final = theta_arr_median
scatter_vec_final = scatter_vec_median




; Write it all
med_outstr = {lambda_mask:lambda_mask_final, $
              scatter_vec:scatter_vec_final, $
              theta_arr:theta_arr_final}
    
mwrfits, med_outstr, median_outfile, /create




;skip_training:

;;; Now, the test step
; Using the theta and scatter, we optimize the labels for each
; spectrum.

;if skip_opt then begin
    training_str = mrdfits(median_outfile,1)
    lambda_mask = training_str.lambda_mask
    scatter_vec = training_str.scatter_vec
    theta_arr = training_str.theta_arr
;endif

; data_cross = spectra_ol[*,cross_idx]
; err_data_cross = error_ol[*,cross_idx]
; mask_data_cross = mask_ol[*,cross_idx]
;data_cross = spectra_ol
data_cross = smooth_spec_ol
err_data_cross = smooth_err
mask_data_cross = mask_ol




;FIXME (this is true, but inflexible)
ncross = nspec_total
cross_idx = lindgen(nspec_total)
teff_train = teff_all
feh_train = feh_all
logg_train = logg_all
vsini_train = vsini_all

;;; Begin loop through spectra
; Initialize

common test, test_matrix, test_sigma, test_data, model_type

model_type = model_names[where(model_settings eq 1)]

nlabels_test = nlabels_set

test_status_vec = bytarr(ncross)
test_label_results = dblarr(nlabels_test, ncross) 
test_label_errors = dblarr(nlabels_test, ncross) 

for snum = 0, ncross-1 do begin

    ; Grab the spectrum
    test_data = data_cross[*,snum]

    test_err = err_data_cross[*,snum]
    
    ; Mask wavelengths that are masked in the spectrum, and those that
    ; failed in training
    test_mask = mask_data_cross[*,snum]
    test_mask_bin = (((mask_dec and test_mask) + lambda_mask) < 1)
    test_mask_idx = where(test_mask_bin ne 0, nmask_test, comp=test_unmask_idx, $
                          ncomp=nunmask_test)
    
    ; Jack up error for masked pixels
    test_err_masked = test_err
    ; if nmask_test gt 0 then begin
;         test_err_masked[test_mask_idx] = 1d8
;     endif

    ; Combine errors
    test_sigma = sqrt(test_err_masked^2 + scatter_vec)


    ; Test_matrix (design matrix) is theta_arr
    test_matrix = theta_arr
    
    ; Make covariance matrix for returning label error estimates
    ; It's [ M^T C^-1 M ]^-1
    
    inverse_covar = diag_matrix(1d0/(test_sigma^2d0))
        
    dm_trans = transpose(test_matrix)
    
    ; Do the thing in the brackets
    invco_dm = inverse_covar ## test_matrix
    thing_in_brackets = dm_trans ## invco_dm
    
    ; Get the output covariance matrix
    covar_matrix = invert(thing_in_brackets, covar_stat, /double)

    ;;; We need to use amoeba to optimize the labels
    ; INPUTS: theta_arr, test_sigma, test_data, guess (a vector of
    ; nlabels_set)
    ; FIXME: Add covariance matrix to common block (test) in order to
    ; get label errors back from the test step


    ftol_test = 1d-5
    ;test_guess_all = [mean(vsini_train), mean(teff_train), mean(feh_train), $
                      ;mean(logg_train)]

    test_guess_all = replicate(0d0, n_elements(label_set_idx))
    ;test_scale_all = [1000d0, 1d0, 2d0, 10d0]
    test_scale_all = replicate(1d0, n_elements(label_set_idx))

    test_guess = test_guess_all[label_set_idx]
    test_scale = test_scale_all[label_set_idx]
    nmax_test = 5000

    ; Trying no AMOEBA
    test_constant_vec = reform(theta_arr[0,*])

    test_matrix_rgs = theta_arr[1:*,*]

    ; Remove the constant offset from test_data
    test_yvals = test_data - test_constant_vec
    
    test_results = regress(test_matrix_rgs, test_yvals, measure_errors=test_sigma, $
                             sigma=label_errors, status=test_status)

    test_results = reform(test_results)
;     test_results = amoeba3(ftol_test, function_name='test_step', $
;                              p0=test_guess, scale=test_scale, $
;                              function_value=fval_test, $
;                              nmax=nmax_test, ncalls=ncalls_test)
    


    

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


; This is pretty dumb
;Vsini
if label_settings[0] eq 1 then begin

    vsini_label_idx = where(label_set_idx eq 0, noops)
    if noops eq 0 then message, "Problem with vsini idx"

    test_vsini = (reform(test_label_results[vsini_label_idx,*])*2d0*vsini_sig) + vsini_mean
    test_e_vsini = (reform(test_label_errors[vsini_label_idx,*])*2d0*vsini_sig)
    res_vsini = test_vsini - vsini_ol
    vsini_true = 1
endif else vsini_true = 0

; Teff
if label_settings[1] eq 1 then begin

    teff_label_idx = where(label_set_idx eq 1, noops)
    if noops eq 0 then message, "Problem with teff idx"

    test_teff = (reform(test_label_results[teff_label_idx,*])*2d0*teff_sig) + teff_mean
    test_e_teff = (reform(test_label_errors[teff_label_idx,*])*2d0*teff_sig)
    res_teff = test_teff - teff_ol
    teff_true = 1
endif else teff_true = 0

; Feh
if label_settings[2] eq 1 then begin

    feh_label_idx = where(label_set_idx eq 2, noops)
    if noops eq 0 then message, "Problem with feh idx"

    test_feh = (reform(test_label_results[feh_label_idx,*])*2d0*feh_sig) + feh_mean
    test_e_feh = (reform(test_label_errors[feh_label_idx,*])*2d0*feh_sig)
    res_feh = test_feh - feh_ol
    feh_true = 1
endif else feh_true = 0

; Logg
if label_settings[3] eq 1 then begin

    logg_label_idx = where(label_set_idx eq 3, noops)
    if noops eq 0 then message, "Problem with logg idx"

    test_logg = (reform(test_label_results[logg_label_idx,*])*2d0*logg_sig) + logg_mean
    test_e_logg = (reform(test_label_errors[logg_label_idx,*])*2d0*logg_sig)
    res_logg = test_logg - logg_ol
    logg_true = 1
endif else logg_true = 0

; NON DET CONTAMINATION
if label_settings[0] eq 1 then begin

    nd_true_idx = where(vsini_ol le 5, nndtrue)
    nd_fit_true_idx = where(test_vsini[nd_true_idx] le 5, nndfittrue, comp=nd_fit_false_idx, ncomp=nndfitfalse)

    nsuccess = double(nndfittrue)/double(nndtrue)

    contamination = double(nndfitfalse)/(double(nndfittrue)+double(nndfitfalse))
endif



; Mark failed guys in red
failed_idx = where(test_status_vec ne 0, nfail_test)


;window, 1, tit="RESULTS!!!", xs=1500, ys=900
plotdir = '/home/stgilhool/Vsini_ML/cannon_test/plots_cannon_test/'
plotfile = plotdir + 'cannon_test_'+test_number+'.eps'

psopen, plotfile, /encaps, /color, xs=10, ys=8, /inches
setcolors, /system_var, /color, /silent
!p.color = !black

!p.multi = [0,nlabels_set, nlabels_set]


;window, 0, xs = 1600, ys=900

;xyouts, 0.5, 0.95, model_type, alignment=1, charsize=3
;xyouts, 0.05, 0.5, "CANNON - CKS", orientation=90, alignment=1, charsize=3



nlabels_tot = n_elements(label_names) 

for lab_i = 0, nlabels_tot - 1 do begin

    
    if label_settings[lab_i] eq 1 then begin
        ; Make color bar for lab_i, if it's included
        case lab_i of

            0: begin
                mag_vec = vsini_ol
                cbar_title = "Vsini (CKS)"
            end

            1: begin
                mag_vec = teff_ol
                cbar_title = "Teff (CKS)"
            end

            2: begin
                mag_vec = feh_ol
                cbar_title = "[Fe/H] (CKS)"
            end

            3: begin
                mag_vec = logg_ol
                cbar_title = "Logg (CKS)"
            end

            
            else: message, "Error in label settings jawn"
                
        endcase
        

    endif else continue

    ; VSINI
    if vsini_true then begin
        ; vsini_plot = scatterplot(vsini_cross, res_vsini, symbol='circle', $
;                                 sym_size=0.2, xtit="Vsini (CKS)", $
;                                 ytit=textoidl('\delta Vsini'), $
;                                  magnitude=mag_vec)
        
        ;         rline_vsini = plot(vsini_plot.xrange, [0,0],
        ;         /overplot, /thick)
        sg_scatter_plot, vsini_ol, res_vsini, test_e_vsini, $
          colors=bytscl(mag_vec, top=255-12), xtit=textoidl('vsini (CKS)'), $
          ytit=textoidl('\Delta vsini')
        oplot, minmax(vsini_ol), replicate(0d0, 2), linest=2, /thick
        
    endif



; TEFF
    if teff_true then begin
        ; teff_plot = scatterplot(teff_cross, res_teff, symbol='circle', $
;                                 sym_size=0.2, 

;                                 magnitude=mag_vec)
        
;         rline_teff = plot(teff_plot.xrange, [0,0], /overplot, /thick)
        sg_scatter_plot, teff_ol, res_teff, sqrt(test_e_teff^2), $
          colors=bytscl(mag_vec, top=255-12), xtit=textoidl('T_{eff} (CKS)'), $
          ytit=textoidl('\Delta T_{eff}')
        oplot, minmax(teff_ol), replicate(0d0, 2), linest=2, /thick


    endif

    ; FEH
    if feh_true then begin
        ; feh_plot = scatterplot(feh_cross, res_feh, symbol='circle', $
;                                sym_size=0.2, xtit="[Fe/H] (CKS)", $
;                                ytit=textoidl('\delta [Fe/H]'), $
;                                magnitude=mag_vec)
        
;         rline_feh = plot(feh_plot.xrange, [0,0], /overplot, /thick)
        sg_scatter_plot, feh_ol, res_feh, test_e_feh, $
          colors=bytscl(mag_vec, top=255-12), xtit=textoidl('[Fe/H] (CKS)'), $
          ytit=textoidl('\Delta [Fe/H]')
        oplot, minmax(feh_ol), replicate(0d0, 2), linest=2, /thick
            
    endif

    ; LOGG
    if logg_true then begin
        ; logg_plot = scatterplot(logg_cross, res_logg, symbol='circle', $
;                                 sym_size=0.2, xtit="Logg (CKS)", $
;                                 ytit=textoidl('\delta logg'), $
;                                 magnitude=mag_vec)
        
        ;         rline_logg = plot(logg_plot.xrange, [0,0],
        ;         /overplot, /thick)
        sg_scatter_plot, logg_ol, res_logg, test_e_logg, $
          colors=bytscl(mag_vec, top=255-12), xtit=textoidl('logg (CKS)'), $
          ytit=textoidl('\Delta logg')
        oplot, minmax(logg_ol), replicate(0d0, 2), linest=2, /thick

    endif

    
;    c = colorbar(orientation=1, title=cbar_title)

endfor

if label_settings[0] eq 1 then begin
    xyouts, 0.5, 0.95, "Vsini Non-detection Contamination = "+strtrim(sigfig(contamination,3),2),/alignment, /normal
endif

psclose
setcolors, /system_var, /color, /silent
!p.background=!black


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

stop

end
