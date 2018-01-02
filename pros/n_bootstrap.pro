;+
; NAME: n_bootstrap
;
;
;
; PURPOSE: Return the number of stars in the full APOGEE sample
; matching the input conditions
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE: n = n_bootstrap(c0,[c1,c2,c3,...c9],
;                                   ALLSTAR=allstar, RET_IDX=ret_idx)
;
;
;
; INPUTS:
;  c0: conditional statement string, referencing a field in allstar
;  (ie. allstar.teff gt 5000)
;  

;
; OPTIONAL INPUTS:
;  c1-9: more conditionals (in same format as c0). Will be joined with
;  an AND statement
;
;
; KEYWORD PARAMETERS:
;  ALLSTAR: structure containing the APOGEE allStar file. Will be
;  read-in, if not supplied, but takes 10's of seconds
;
;
;
; OUTPUTS:
;  N: number of stars that fit the conditionals
;
;
; OPTIONAL OUTPUTS:
;  RET_IDX: vector of indices to the stars that fit the conditionals
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

function n_bootstrap, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, ALLSTAR=allstar, RET_IDX=ret_idx

;;; Check input 
if n_elements(ret_idx) eq 0 then ret_idx = []
; (must be no more than 10 string conditions)
n_cond = n_params() 

if ((n_cond eq 0) or (n_cond gt 10)) then message, "You must enter 1-10 conditional statments"

; (readin allstar file, if not passed as input)
if n_elements(allstar) eq 0 then begin
    asfile = '/home/stgilhool/APOGEE/APOGEE_data/allStar-l31c.2.fits'
    allstar = mrdfits(asfile, 1)
endif

;;; Manipulate input variables into a conditional statement to be
;;; passed to where()
; Make string array of variable names
cond_var_arr = 'c'+strtrim(indgen(n_cond),2)

; Loop through to create the string conditional
conditional_arr = []

foreach cond_var, cond_var_arr do begin

    cond_string_i = scope_varfetch(cond_var, level=0)
    conditional_arr = [conditional_arr, cond_string_i]

endforeach

; Join the strings into one
if n_cond eq 1 then conditional_statement = conditional_arr[0] $
  else conditional_statement = strjoin(conditional_arr, ' and ')

;;; Find the stars that meet the condtions
void = execute('return_index = where('+conditional_statement+',n_match)')

if n_match eq 0 then message, "No stars meet the given conditions" $
  else ret_idx = return_index

return, n_match
    
end
