;;; Procedure to set and output lasso_select settings

pro lasso_settings, runname, WRITE=write, READ=read

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
    endcase

endif else if rundir_check ne 1 then message, "Some weird error with run directory"

; Now, check on the settings file, 'runfile'
runfile_name = runname+'.fits'
runfile = rundir+runfile_name

runfile_chk = file_test(runfile)

; If user doesn't specify READ and WRITE, figure out what to do
if n_elements(write) eq 0 and n_elements(read) eq 0 then begin
; if the file exists already, then READ
if runfile_chk eq 1 then begin
    readset = 1 
    writeset = 0
; if the file doesn't exist already, then WRITE
endif else if runfile_chk eq 0 then begin
    readset = 0
    writeset = 1
endif else message, "Some weird error with runfile_chk"

; If both are specified, set both
endif else if n_elements(write) ne 0 and n_elements(read) ne 0 then begin
    writeset = 1
    readset = 1
; Otherwise, zero the not-set one
endif else begin
    writeset = keyword_set(write)
    readset = keyword_set(read)
    if writeset + readset ne 1 then message, "Whoops: check use of keyword_set!"
endelse



end
    

