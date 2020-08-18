;***************************************
; PROGRAM DESCRIPTION: Create timeseries of CESM ice volume
; INPUT DATA: CESM timeseries
; OUTPUT DATA: plot
; CREATOR: Alice DuVivier - Jan 2019
;***************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;***************************************
begin
;;;;;;;;;;
; Manual inputs or from cshell script
;;;;;;;;;;
; ******************************************************
varcode = "albedo"
hemi = "NH"     ; choose hemisphere - "NH" or "SH"
; ******************************************************
; BEGIN SCRIPT

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
;data_dir = "./"
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_ice/analysis/version_9/albedo_annual_cycle/"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; load albedo difference from atmosphere
; over: ice >= 15% concentration + open water
f1 = addfile(data_dir+"atmosphere_albedo_annual_cycle/CMIP6_b_bw_"+hemi+"_piControl_albedo_diff_annual_cycle_seaice_only_for_plot.nc","r")
diff_atm = f1->diff
pcnt_diff_atm = f1->pcnt_diff
delete(f1)

; load albedo difference from ice
; over: ice >= 15% concentration only 
; (DOES NOT INCLUDE IMPACT FROM ICE FRACTION CHANGE)
f2 = addfile(data_dir+"ice_albedo_annual_cycle/CMIP6_b_bw_"+hemi+"_piControl_albedo_diff_annual_cycle_for_plot.nc","r")
diff_ice = f2->diff
pcnt_diff_ice = f2->pcnt_diff
delete(f2)

;;;;;;;;;;
; Combine for plotting
;;;;;;;;;;
nmons = 12
diff_all = new((/4,nmons+1/),float)
diff_all(0:1,0:nmons-1) = diff_ice
diff_all(2:3,0:nmons-1) = diff_atm
diff_all(0:1,nmons) = diff_ice(:,0)
diff_all(2:3,nmons) = diff_atm(:,0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_diff_annual_cycle_both_ice_atm"

xdata = ispan(1, nmons+1, 1)
label_inds = xdata
label_str = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)

if (hemi .eq. "NH")
   pref = "a) Albedo Difference"
   ymin = -0.05
   ymax = 0.01
   yminor = 3
   x0 = 6.0
   x1 = 4.0
   x2 = x1
   y0 = 0.005
   y1 = -0.035
   y2 = -0.04 
end if
if (hemi .eq. "SH")
   pref = "d) Antarctic Difference~C~(CAM6-WACCM6)"
   ymin = -0.005
   ymax = 0.001
   yminor = 1
   x0 = 3.0
   y0 = 0.0005
end if

units = ""

;;;;;;;;;;;;;
; XY Cross Section Plot Resources
;;;;;;;;;;;;;
; set resources used for all plots
  opts                             = True              ; plot mods desired
  opts@gsnFrame                    = False             ; Do not draw plot 
  opts@gsnDraw                     = False             ; Do not advance frame
  opts@gsnMaximize                 = False             ; maximize plot
  opts@vpWidthF                    = 1.25  ;1.75               ; set aspect ratio
  opts@vpHeightF                   = 0.75 ;0.55

; set title resources
  opts@tiMainOn                    = True
  opts@tiMainFontHeightF           = 0.03
  opts@tiMainFontThicknessF        = 1.5
  opts@tmXBOn                      = True
  opts@tmXBMode                    = "Explicit"
  opts@tmXBValues                  = label_inds
  opts@tmXBLabels                  = label_str
  opts@trXMinF                     = 1
  opts@trXMaxF                     = nmons+1
  opts@tiXAxisString               = "Month"
  opts@tiYAxisString               = units
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@tmYLMinorPerMajor           = yminor
  opts@tmYRMinorPerMajor           = yminor
  opts@trYMinF                     = ymin
  opts@trYMaxF                     = ymax
  opts@gsnYRefLine                  = 0.0

; line resources
  opts@tiMainOn = False
  opts@xyMarkLineMode = "MarkLines"
  opts@xyLineColors = (/"black", "black", "grey40", "grey40"/)
  opts@xyLineThicknessF = 3.5
  opts@xyDashPatterns = (/0,0,2,2/)
  opts@xyMarkers = (/4,16,4,16/)
  opts@xyMarkerColors = opts@xyLineColors

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2
  txres_t = txres
  txres_t@txFontColor = "black"
  txres_b = txres
  txres_b@txFontColor = "black"
  txres_bw = txres
  txres_bw@txFontColor = "grey40"

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 3                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 3                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True	       ; fill the page

;;;;;;;;;;;;;
; Make plots
;;;;;;;;;;;;;
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)    

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,diff_all,opts)
   dum0 = gsn_add_text(wks,lplot,pref,x0,y0,txres_t)
   dum1 = gsn_add_text(wks,lplot,"Albedo (ice only)",x1,y1,txres_b)
   dum2 = gsn_add_text(wks,lplot,"Albedo (ice + open water)",x2,y2,txres_bw)
   plot(0) = lplot

;;;;;;;;;;;;;
; Final Plots
;;;;;;;;;;;;;
  ; make panel plot
  ;pnlres@txString  = hemi+" Monthly Sea Ice Extent (mean over entire PIctrl)"
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Figure done, yay!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
