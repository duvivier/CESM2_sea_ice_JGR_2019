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
varcode = "sithick"
hemi = "NH"
; ******************************************************
; BEGIN SCRIPT
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "./"

; historical period: 1979-2014

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Load data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; load historical
f1 = addfile(data_dir+"CMIP6_b_bw_"+hemi+"_hist_fractions_for_histogram.nc","r")
; load bins
levels_plot_hist = f1->levels_plot
; load data from historical file
b_spring_area_frac_hist = f1->b_spring_area_frac
bw_spring_area_frac_hist = f1->bw_spring_area_frac
b_fall_area_frac_hist = f1->b_fall_area_frac
bw_fall_area_frac_hist = f1->bw_fall_area_frac
delete(f1)

; load pi control
f1 = addfile(data_dir+"CMIP6_b_bw_"+hemi+"_pi_fractions_for_histogram.nc","r")
; load bins
levels_plot_pi = f1->levels_plot
; load data from historical file
b_spring_area_frac_pi = f1->b_spring_area_frac
bw_spring_area_frac_pi = f1->bw_spring_area_frac
b_fall_area_frac_pi = f1->b_fall_area_frac
bw_fall_area_frac_pi = f1->bw_fall_area_frac
delete(f1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting PDF of "+varcode)

levels_plot = levels_plot_hist
label_inds = ispan(0,6,1)
label_str = flt2string(int2flt(label_inds))
if (hemi .eq. "NH")
  seas_str1 = "a) Arctic Winter - JFM"
  seas_str2 = "b) Arctic Summer - JAS"
end if
if (hemi .eq. "SH")
  seas_str1 = "c) Antarctic Summer - JFM"
  seas_str2 = "d) Antarctic Winter - JAS"
end if

; locations for text
x0 = 2.2
x1 = 4.6
x2 = 4.3
x3 = 5.15
x4 = 4.85
y0 = 0.38
y1 = 0.35
y2 = 0.32
y3 = 0.29
y4 = 0.26

; set outdir and fout names
   outdir = "./"

; set general resources
  res                     = True          ; plot mods desired
  res@gsnDraw             = False         ; don't draw the plot yet
  res@gsnFrame            = False         ; or the frame

; bar chart resources
  res@gsnXYBarChart       = True            ; Create bar plot
  res@gsnXYBarChartOutlineOnly = True 
  res@xyLineThicknessF    = 5.
  res@xyDashPattern       = 0
  res@xyMonoLineColor      = False

; labeling resources
  res@tmXBMode            = "Explicit"         ; explicit labels
  res@tmXBValues          = label_inds
  res@tmXBLabels          = label_str
  res@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  res@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  res@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  res@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  res@trXMinF             = 0
  res@trXMaxF             = 6.0
  res@trYMinF             = 0
  res@trYMaxF             = 0.4
  res@tiYAxisString       = "Fraction of total ice area"
  res@tiXAxisString       = "ice thickness bins [m]"

; text resources
  txres = True
  txres@txFontHeightF = 0.025 ;0.03
  txres@txFontThicknessF = 2

; CAM resources - PI
  res_b_pi = res
  res_b_pi@xyLineColors = (/"grey50"/)
  txres_b_pi = txres
  txres_b_pi@txFontColor = "grey50"

; WACCM resources- PI
  res_bw_pi = res
  res_bw_pi@xyLineColors = (/"lightpink3"/)
  txres_bw_pi = txres
  txres_bw_pi@txFontColor = "lightpink3"

; CAM resources - Hist
  res_b_h = res
  res_b_h@xyLineColors = (/"black"/)
  txres_b_h = txres
  txres_b_h@txFontColor = "black"

; WACCM resources- Hist
  res_bw_h = res
  res_bw_h@xyLineColors = (/"firebrick"/)
  txres_bw_h = txres
  txres_bw_h@txFontColor = "firebrick"

;; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@lbOrientation              = "vertical"
   pnlres@gsnPanelFigureStringsBackgroundFillColor = "white"  ; makes background of label box white
  pnlres@gsnPanelFigureStringsFontHeightF = 0.015   ; text height of label
  pnlres@amJust                     = "TopRight"                 ; put label in top right corner
  pnlres@gsnPanelFigureStringsPerimOn = False ; no box around labels

;;;;;;;;;;;;;
; Make regular plots
;;;;;;;;;;;;;
print("Plotting seasonal histograms")
; set info for this plot
   fout = outdir+"CMIP6_b_bw_"+hemi+"_"+varcode+"_histogram"

  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,"wh-bl-gr-ye-re") 

; do individual plots
plot = new(2,graphic)
s = 0
bplot = gsn_csm_xy(wks,levels_plot,b_spring_area_frac_hist,res_b_h)
oplot1 = gsn_csm_xy(wks,levels_plot,bw_spring_area_frac_hist,res_bw_h)
overlay(bplot,oplot1)
oplot2 = gsn_csm_xy(wks,levels_plot,b_spring_area_frac_pi,res_b_pi)
overlay(bplot,oplot2)
oplot3 = gsn_csm_xy(wks,levels_plot,bw_spring_area_frac_pi,res_bw_pi)
overlay(bplot,oplot3)
dum0 = gsn_add_text(wks,bplot,seas_str1,x0,y0,txres_b_h)
dum1 = gsn_add_text(wks,bplot,"CAM6 - Historical",x1,y1,txres_b_h)
dum2 = gsn_add_text(wks,bplot,"WACCM6 - Historical",x2,y2,txres_bw_h)
dum3 = gsn_add_text(wks,bplot,"CAM6 - PI",x3,y3,txres_b_pi)
dum4 = gsn_add_text(wks,bplot,"WACCM6 - PI",x4,y4,txres_bw_pi)
plot(s) = bplot

s = 1
bplot = gsn_csm_xy(wks,levels_plot,b_fall_area_frac_hist,res_b_h)
oplot1 = gsn_csm_xy(wks,levels_plot,bw_fall_area_frac_hist,res_bw_h)
overlay(bplot,oplot1)
oplot2 = gsn_csm_xy(wks,levels_plot,b_fall_area_frac_pi,res_b_pi)
overlay(bplot,oplot2)
oplot3 = gsn_csm_xy(wks,levels_plot,bw_fall_area_frac_pi,res_bw_pi)
overlay(bplot,oplot3)
;dum1 = gsn_add_text(wks,bplot,"CAM6 - Historical",x1,y1,txres_b_h)
;dum2 = gsn_add_text(wks,bplot,"WACCM6 - Historical",x2,y2,txres_bw_h)
;dum3 = gsn_add_text(wks,bplot,"CAM6 - PI",x3,y3,txres_b_pi)
;dum4 = gsn_add_text(wks,bplot,"WACCM6 - PI",x4,y4,txres_bw_pi)
dum5 = gsn_add_text(wks,bplot,seas_str2,x0,y0,txres_b_h)
plot(s) = bplot

; add text overlay

; panel plot
print("Paneling")
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,2/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Done with all plots. Hooray!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
