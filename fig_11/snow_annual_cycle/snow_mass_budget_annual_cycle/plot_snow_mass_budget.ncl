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
varcode = "net_mass_budget"
hemi = "NH"     ; choose hemisphere - "NH" or "SH"
; ******************************************************
; BEGIN SCRIPT

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_ice/analysis/version_9/snow_annual_cycle/snow_mass_budget_annual_cycle/"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f1 = addfile(data_dir+"/CMIP6_b_bw_"+hemi+"_piControl_snow_mass_budget_terms_for_plot.nc","r")
diff_plots = f1->diff_plots
b_data_all = f1->b_data_all
bw_data_all = f1->bw_data_all
delete(f1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"

xdata = ispan(1, 12, 1)
label_inds = xdata
label_str = flt2string(int2flt(xdata))

units = "cm/day"

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
  opts@trXMaxF                     = 12
  opts@tiXAxisString               = "Month"
  opts@tiYAxisString               = units
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@gsnYRefLine                 = 0.0

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2

; resources for b
  opts_b = opts
  opts_b@tiMainOn = False
  opts_b@xyMarkLineMode = "MarkLines"
  opts_b@xyLineColors = (/"black","blue","firebrick","darkorange","black","blue","firebrick","darkorange"/)
  opts_b@xyMarkerColors = opts_b@xyLineColors
  opts_b@xyMarkers = (/1,1,1,1,16,16,16,16/)
  opts_b@xyLineThicknesses = (/5.,4.,4.,4.,1.,1.,1.,1./)
  opts_b@xyDashPattern = 0

; text
  title1 = "Net"
  txres_1 = txres
  txres_1@txFontColor = "black"
  title2 = "Snowfall"
  txres_2 = txres
  txres_2@txFontColor = "blue"
  title3 = "Melt"
  txres_3 = txres
  txres_3@txFontColor = "firebrick"
  title4 = "Evap/Subl"
  txres_4 = txres
  txres_4@txFontColor = "darkorange"

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
; Info for labels
;;;;;;;;;;;;;
if (hemi .eq. "NH")
   pref = "Snow Mass Budget"
   ymin = -1.0
   ymax = 0.5
   yminor = 2
   x0 = 6.0
   x1 = 3.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 0.4
   y1 = -0.3
   y2 = -0.4
   y3 = -0.5
   y4 = -0.6
end if
if (hemi .eq. "SH")
   pref = "Snow Mass Budget"
   ymin = -0.8
   ymax = 0.5
   yminor = 3
   x0 = 6.0
   x1 = 4.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 0.45
   y1 = -0.3
   y2 = -0.4
   y3 = -0.5
   y4 = -0.6
end if

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor

;;;;;;;;;;;;;
; Make b case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_"+hemi+"_piControl_snow_mass_budget"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic) 

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,b_data_all,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+" CAM6",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   plot(0) = lplot

;;;;;;;;;;;;;
; Final Plots
;;;;;;;;;;;;;
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/plot,dum0,dum1,dum2,dum3,dum4/])
print("Figure done, yay!")

;;;;;;;;;;;;;
; Make bw case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_bw_"+hemi+"_piControl_snow_mass_budget"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic) 

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,bw_data_all,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+" WACCM6",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   plot(0) = lplot

;;;;;;;;;;;;;
; Final Plots
;;;;;;;;;;;;;
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/plot,dum0,dum1,dum2,dum3,dum4/])
print("Figure done, yay!")

;;;;;;;;;;;;;
; Make difference plots
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_snow_mass_budget_diffs"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

if (hemi .eq. "NH")
   pref = "c) Arctic Snow on Ice Mass Budget"
   ymin = -0.04
   ymax = 0.04
   yminor = 4
   x0 = 6.0
   x1 = 10.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 0.035
   y1 = 0.028
   y2 = 0.024
   y3 = 0.020
   y4 = 0.016
end if
if (hemi .eq. "SH")
   pref = "b) Mass Budget Difference"
   ymin = -0.05
   ymax = 0.05
   yminor = 3
   x0 = 3.5
   x1 = 7.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 0.04
   y1 = 0.045
   y2 = 0.039
   y3 = 0.033
   y4 = 0.027
end if

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,diff_plots,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+" (CAM6-WACCM6)",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   plot(0) = lplot

;;;;;;;;;;;;;
; Final Plots
;;;;;;;;;;;;;
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Figure done, yay!")

return
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_snow_mass_budget_diffs"

xdata = ispan(1, 12, 1)
label_inds = xdata
label_str = flt2string(int2flt(xdata))

if (hemi .eq. "NH")
   pref = "c) Arctic Snow Mass Budget"
   ymin = -0.05
   ymax = 0.05
   yminor = 4
   x0 = 6.0
   x1 = 8.0
   x2 = x1
   x3 = 10.0
   x4 = x3
   y0 = 0.023
   y1 = 0.018
   y2 = 0.014
   y3 = 0.018
   y4 = 0.014
end if
if (hemi .eq. "SH")
   pref = "c) Antarctic Mass Budget"
   ymin = -0.025
   ymax = 0.055
   yminor = 3
   x0 = 3.5
   x1 = 7.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 10.0
   x6 = x5
   x7 = x5
   x8 = x5
   x9 = x5
   y0 = 0.04
   y1 = 0.045
   y2 = 0.039
   y3 = 0.033
   y4 = 0.027
   y5 = 0.045
   y6 = 0.039
   y7 = 0.033
   y8 = 0.021
   y9 = 0.027
end if


units = "cm/day"
;units = "kg m~S~-2~N~ s~S~-1~N~"

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
  opts@trXMaxF                     = 12
  opts@tiXAxisString               = "Month"
  opts@tiYAxisString               = units
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@trYMinF                     = ymin
  opts@trYMaxF                     = ymax
  opts@tmYLMinorPerMajor           = yminor
  opts@tmYRMinorPerMajor           = yminor
  opts@gsnYRefLine                 = 0.0

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; text resources
  txres = True
  txres@txFontHeightF = 0.03 ;0.02
  txres@txFontThicknessF = 2

; resources for b
  opts_b = opts
  opts_b@tiMainOn = False
  opts_b@xyMarkLineMode = "MarkLines"
  opts_b@xyLineColors = (/"black","blue","firebrick","darkorange","black","blue","firebrick","darkorange"/)
  opts_b@xyMarkerColors = opts_b@xyLineColors
  opts_b@xyMarkers = (/1,1,1,1,16,16,16,16/)
  opts_b@xyLineThicknesses = (/5.,4.,4.,4.,1.,1.,1.,1./)
  opts_b@xyDashPattern = 0

; text
  title1 = "Net"
  txres_1 = txres
  txres_1@txFontColor = "black"
  title2 = "Snowfall"
  txres_2 = txres
  txres_2@txFontColor = "blue"
  title3 = "Melt"
  txres_3 = txres
  txres_3@txFontColor = "firebrick"
  title4 = "Evap/Subl"
  txres_4 = txres
  txres_4@txFontColor = "darkorange"

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
   lplot = gsn_csm_xy(wks,xdata,diff_plots,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+" (CAM6-WACCM6)",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
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
