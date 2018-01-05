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
lasso_base_path = '/home/stgilhool/Vsini_ML/lasso_select/'
lasso_trials_path = lasso_base_path+'Test_Runs/'

; Find the directory containing the info for this run, if it exists
rundir = lasso_trials_path+runname+'/'

rundir_chk = file_test(rundir, /directory)

; ask to make the directory if it doesn't exist, and exit on error
if rundir_chk eq 0 begin
    
    ; create directory?
    mkdir_opt = ''
    prompt_for_mkdir:
    read, mkdir_opt, $
      PROMPT="Directory for "+runname+" doesn't exist. Create? (y) or (n)"
    
    case mkdir_opt of
        'y' || 'Y': begin
            print, ''
            print, "Making directory: "+rundir
            print, ''
            file_mkdir, rundir
        end
        'n' || 'N': begin
            print, ''
            message, "Not making directory. Exiting..."
        end
        else: begin
            print, ''
            print, "Input not recognized. Please type 'y' or 'n'."
            print, ''
            goto, prompt_for_mkdir
        end
    endcase

endif else if rundir_check ne 1 then message, "Some weird error with run directory"

; Create directory structure, if necessary

; Create default structure

; Update structure with optional inputs, if supplied

; Ask user to review the settings

; Convert structure to JSON object

; Write JSON object

end

