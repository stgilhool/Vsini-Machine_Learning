; Function to partition data into training set, cross-validation set, etc.
function ml_partition_data, param_arr, ntraining, RANDOM=random, SPAN=span, DEBUG=debug

;on_error, 2

; Keywords
if n_elements(random) eq 0 then random = 0 else random = 1
if not keyword_set(span) then span = 0 else span = 1 
if n_elements(debug) eq 0 then debug = 0 else debug = 1

; Arguments
if n_params() ne 2 then begin
    print, "Syntax - ML_PARTITION_DATA, param_arr, ntraining,"
    print, "                            RANDOM=random, SPAN=span, "
    print, "                            DEBUG=debug"
    return, !null
endif

; Initialize
nparam = n_elements(param_arr[0,*]) 
ntotal = n_elements(param_arr[*,0]) 

; Set up partitioning
if n_elements(ntraining) le 0 then message, "ntraining is a required parameter" $
else if n_elements(ntraining) gt 1 then begin
    ; If ntraining is an array, create a structure containing fields
    ; for each set
    nsets = n_elements(ntraining) 
    npts_vec = ntraining

endif else begin
    ; If ntraining is a scalar, it's the number of points to use for
    ; the training set.  We will by default divide the remaining into
    ; a cross-validation set and a test set
    nsets = 3
    
    nremaining = ntotal - ntraining

    ; number of cross-validation points
    nxval = nremaining/2
    ; number of test points
    ntest = nremaining - nxval
    
    npts_vec = [ntraining, nxval, ntest]

endelse

; Construct the structure to be populated later
set_struct = []

; fields are called SET0, SET1, etc.
keys = transpose(transpose(replicate('SET',nsets))+ $
                 transpose(strtrim(lindgen(nsets),2)))

for set_i = 0, nsets-1 do begin
    
    set_struct = create_struct(set_struct, keys[set_i], $
                               replicate(0L, npts_vec[set_i]))

endfor

set_idx_counter = replicate(0L, nsets)

;;;; Populate the sets with points
; There are 3 possible schemes:
; 1) Default - points are distributed so that the training set gets
;              the maximum val point of the first parameter, and the
;              mininum. Each remaining set gets the remaining max and
;              min. The process is repeated with the remaining
;              parameters until all points are assigned to a set
; 2) RANDOM - if random is set, the points are distributed randomly
; 3) SPAN - like default, but only the first max and min are ensured
;           per set.  The remaining points are filled in randomly
; ONLY DEFAULT IS OPERATIONAL NOW


; Loop trough parameter arrays to get the min and max into the
; training set
param_setting_arr = double(param_arr)
terminate = 0

while terminate eq 0 do begin
    ; Go parameter-by-parameter, assigning max and then min to each set
    for parami = 0, nparam - 1 do begin

        for set_i = 0, nsets-1 do begin
        
            par_vec = param_setting_arr[*,parami]
            
            pmin = min(par_vec, pminidx, /nan)
            pminidx = pminidx[0]

            pmax = max(par_vec, pmaxidx, /nan)
            pmaxidx = pmaxidx[0]

            ; Check if there is room left in the given set
            n_in_set = set_idx_counter[set_i]
            
            ; If there is room in the set, add the maximum
            if n_in_set lt npts_vec[set_i] then begin
                ; Add maximum to the given set
                set_struct.(set_i)[set_idx_counter[set_i]] = pmaxidx
                ; increment the index of the set
                set_idx_counter[set_i] = set_idx_counter[set_i] + 1
                ; NaN all values for that point
                param_setting_arr[pmaxidx, *] = !values.d_nan
            endif else continue
            
            ; Check if there is still room left in the given set
            n_in_set = set_idx_counter[set_i]
            ; If there is room in the set, add the MINIMUM
            if n_in_set lt npts_vec[set_i] then begin
                ; Add maximum to the given set
                set_struct.(set_i)[set_idx_counter[set_i]] = pminidx
                ; increment the index of the set
                set_idx_counter[set_i] = set_idx_counter[set_i] + 1
                ; NaN all values for that point
                param_setting_arr[pminidx, *] = !values.d_nan
            endif else continue

        endfor
        
    endfor
    
    if span then begin

        for set_i = 0, nsets-1 do begin
            ; get the number of guys in the set
            n_in_set = set_idx_counter[set_i]
            ; Fill up the rest of the set randomly
            while n_in_set lt npts_vec[set_i] do begin
                ; generate random index
                defsysv, '!RNG', exists=rng_check
                if rng_check eq 0 then begin
                    defsysv, '!RNG', obj_new('RandomNumberGenerator')
                    print, "Defining new system variable"
                endif 
                rnd_num = !rng -> GetRandomNumbers(1)
                random_idx = rnd_num * (ntotal-1)
                random_idx = round(random_idx[0])
                ; check if it's valid
                if total(finite(param_setting_arr[random_idx,*])) eq nparam then begin
                    ; add to the set
                    set_struct.(set_i)[set_idx_counter[set_i]] = random_idx
                    ; increment the index of the set
                    set_idx_counter[set_i] = set_idx_counter[set_i] + 1
                    ; NaN all values for that point
                    param_setting_arr[random_idx, *] = !values.d_nan
                endif
                ; update how many members are in the set
                n_in_set = set_idx_counter[set_i]
            endwhile
        endfor
        
        ; all points should be assigned, so terminate
        terminate = 1

    endif else if total(finite(param_setting_arr)) eq 0 then terminate = 1

endwhile

return, set_struct
        
end

