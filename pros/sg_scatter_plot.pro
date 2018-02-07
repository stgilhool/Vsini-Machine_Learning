;+
; NAME: sg_scatter_plot
;  
;
;
; PURPOSE: Plot a color-coded scatter plot
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE: sg_scatter_plot, x, y, [error], COLORS=colors, _EXTRA=extra
;
;
;
; INPUTS:
; x - x data
; y - y data
; 
; OPTIONAL INPUTS:
; error - error bar size
;
;
; KEYWORD PARAMETERS:
; colors - vector of type byte, containing color of each data point
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
; Since I use setcolors, /system_variables, be sure to scale colors
; like this:
; input_colors = bytscl(magnitude_vector, top = 255-12)
; sg_scatter_plot, x, y, error, colors=input_colors
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


pro sg_scatter_plot, x, y, error, COLORS=colors, _REF_EXTRA=ex

; Return to user if error
on_error, 2

; Save plot symbol
psym_save = !p.psym
; Set default plot symbol to 8
!p.psym = 8

; Reset default plot symbol and exit on error
catch, error_status

if error_status ne 0 then begin
    !p.psym = psym_save
    help, /last_message, output=msg
    print, msg
    catch, /cancel
    return
endif

; Check number of parameters, and if error is an input
if n_params() eq 2 then begin
    errset=0
endif else if n_params() eq 3 then begin
    errset=1
endif else message, "sg_scatter_plot, x, y [, error], COLORS=colors, _EXTRA=ex"

; If no colors are passed, just plot all same color. If one color
; passed, then make it a vector so function don't fail
if n_elements(colors) eq 0 then colors = replicate(0B, n_elements(x)) else $
  if n_elements(colors) eq 1 then colors = replicate(colors, n_elements(x))

; Match PSYM and SYMSIZE kwds
if n_elements(ex) ne 0 then begin
    exps = stregex(ex, '^ps(y(m)?)?$', /extract, /fold_case)
    exsymsize = stregex(ex, '^syms(i(z(e)?)?)?$', /extract, /fold_case)
endif else begin
    exps = ''
    exsymsize = ''
endelse

; Make Initial Plot
plot, x, y, _extra=ex, /nodata

; Loop through each color and plot the points, with errors, if necessary
foreach color_i, colors, co_idx do begin

    ; Plot colored data
    oplot, [x[co_idx]], [y[co_idx]], color=color_i, _extra=ex
    ; Overplot colored error, if necessary
    if errset then begin
        oploterror, [x[co_idx]], [y[co_idx]], [error[co_idx]], $
          color=color_i, errcolor=color_i, _extra=[exps, exsymsize]
    endif
    
endforeach

; Reset default plot symbol
!p.psym = psym_save

end
