; Pro to generate lasso settings JSON file
pro lasso_config, runname, _EXTRA=ex

; Make sure runname is specified
if n_params() eq 0 then message, "Must specify string <runname>"
; and it's a string
if size(runname, /type) ne 7 then message, "<runname> must be string type"

; Convert any illegal characters in runname
rname = idl_validname(runname, /convert_all)
; Tell user if conversion happened
if strtrim(rname,2) ne strtrim(runname,2) $
  then print, runname, rname, FORMAT='("Converting",A0,"to",A0)'

; Check if runname points to an existing file or not

; Check if runname points to an existing directory or not

; Create directory structure, if necessary

; Create default structure

; Update structure with optional inputs, if supplied

; Ask user to review the settings

; Convert structure to JSON object

; Write JSON object

end

