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
; cesm directory
;data_dir = "./"
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_ice/analysis/version_10/ice_mass_budget_annual_cycle/"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; load ice mass budget
f1 = addfile(data_dir+"/CMIP6_b_bw_"+hemi+"_piControl_mass_budget_terms_for_plot.nc","r")
b_data_all = f1->b_data_all
bw_data_all = f1->bw_data_all
diff_plots = f1->diff_plots
delete(f1)

nmons = 12
; wrap b_data_all and bw_data_all
b_data_all_wrap = new((/9,nmons+1/),float)
b_data_all_wrap(:,0:nmons-1) = b_data_all(:,:)
b_data_all_wrap(:,nmons) = b_data_all(:,0)
bw_data_all_wrap = new((/9,nmons+1/),float)
bw_data_all_wrap(:,0:nmons-1) = bw_data_all(:,:)
bw_data_all_wrap(:,nmons) = bw_data_all(:,0)
diff_plots_wrap = new((/18,nmons+1/),float)
diff_plots_wrap(:,0:nmons-1) = diff_plots(:,:)
diff_plots_wrap(:,nmons) = diff_plots(:,0)


; subset for just biggest terms
diff_plots_sub = new((/10,nmons+1/),float)
diff_plots_sub(0,0:nmons-1) = diff_plots(0,:)   ; net difference
diff_plots_sub(1,0:nmons-1) = diff_plots(1,:)   ; frazil
diff_plots_sub(2,0:nmons-1) = diff_plots(2,:)   ; congelation
diff_plots_sub(3,0:nmons-1) = diff_plots(4,:)   ; top melt
diff_plots_sub(4,0:nmons-1) = diff_plots(5,:)   ; bottom melt
diff_plots_sub(5,0:nmons-1) = diff_plots(9,:)
diff_plots_sub(6,0:nmons-1) = diff_plots(10,:)
diff_plots_sub(7,0:nmons-1) = diff_plots(11,:)
diff_plots_sub(8,0:nmons-1) = diff_plots(13,:)
diff_plots_sub(9,0:nmons-1) = diff_plots(14,:)

; wrap back to january
diff_plots_sub(0,nmons) = diff_plots(0,0)   ; net difference
diff_plots_sub(1,nmons) = diff_plots(1,0)   ; frazil
diff_plots_sub(2,nmons) = diff_plots(2,0)   ; congelation
diff_plots_sub(3,nmons) = diff_plots(4,0)   ; top melt
diff_plots_sub(4,nmons) = diff_plots(5,0)   ; bottom melt
diff_plots_sub(5,nmons) = diff_plots(9,0)
diff_plots_sub(6,nmons) = diff_plots(10,0)
diff_plots_sub(7,nmons) = diff_plots(11,0)
diff_plots_sub(8,nmons) = diff_plots(13,0)
diff_plots_sub(9,nmons) = diff_plots(14,0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"

xdata = ispan(1, 13, 1)
label_inds = xdata
label_str = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)

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
  opts@trXMaxF                     = 13
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
  opts_b@xyLineColors = (/"black","blue","deepskyblue","purple","firebrick","palevioletred1","darkorange","green3","gold","black","blue","deepskyblue","purple","firebrick","palevioletred1","darkorange","green3","gold"/)
  opts_b@xyMarkerColors = opts_b@xyLineColors
  opts_b@xyMarkers = (/1,1,1,1,1,1,1,1,1,16,16,16,16,16,16,16,16,16/)
  opts_b@xyLineThicknesses = (/5.,4.,4.,4.,4.,4.,4.,4.,4.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
  opts_b@xyDashPattern = 0

; text
  title1 = "Net"
  txres_1 = txres
  txres_1@txFontColor = "black"
  title2 = "Frazil Growth"
  txres_2 = txres
  txres_2@txFontColor = "blue"
  title3 = "Congelation Growth"
  txres_3 = txres
  txres_3@txFontColor = "deepskyblue"
  title4 = "Snowice Growth"
  txres_4 = txres
  txres_4@txFontColor = "purple"
  title5 = "Top Melt"
  txres_5 = txres
  txres_5@txFontColor = "firebrick"
  title6 = "Bottom Melt"
  txres_6 = txres
  txres_6@txFontColor = "palevioletred1"
  title7 = "Lateral Melt"
  txres_7 = txres
  txres_7@txFontColor = "darkorange"
  title8 = "Dynamics"
  txres_8 = txres
  txres_8@txFontColor = "green3"
  title9 = "Evap/Subl"
  txres_9 = txres
  txres_9@txFontColor = "gold"

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
   pref = "c) Sea Ice Mass Budget"
   ymin = -2.0
   ymax = 1.0
   yminor = 4
   x0 = 6.5
   x1 = 3.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 10.
   x6 = x5
   x7 = x5
   x8 = x1
   x9 = x5
   y0 = 0.8
   y1 = -0.7
   y2 = -0.9
   y3 = -1.1
   y4 = -1.3
   y5 = -0.7
   y6 = -0.9
   y7 = -1.1
   y8 = -1.5
   y9 = -1.3
end if
if (hemi .eq. "SH")
   pref = "a) Mass Budget"
   ymin = -0.8
   ymax = 0.5
   yminor = 3
   x0 = 6.0
   x1 = 4.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 7.
   x6 = x5
   x7 = x5
   x8 = x1
   x9 = x5
   y0 = 0.45
   y1 = -0.3
   y2 = -0.4
   y3 = -0.5
   y4 = -0.6
   y5 = -0.4
   y6 = -0.5
   y7 = -0.6
   y8 = -0.7
   y9 = -0.7
end if

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor

;;;;;;;;;;;;;
; Make b case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_"+hemi+"_piControl_mass_budget"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic) 

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,b_data_all_wrap,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+" CAM6",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   dum5 = gsn_add_text(wks,lplot,title5,x5,y5,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x6,y6,txres_6)
   dum7 = gsn_add_text(wks,lplot,title7,x7,y7,txres_7)
   dum8 = gsn_add_text(wks,lplot,title8,x8,y8,txres_8)
   dum9 = gsn_add_text(wks,lplot,title9,x9,y9,txres_9)
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

delete([/plot,dum0,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9/])
print("Figure done, yay!")

;;;;;;;;;;;;;
; Make bw case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_bw_"+hemi+"_piControl_mass_budget"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic) 

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,bw_data_all_wrap,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+" WACCM6",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   dum5 = gsn_add_text(wks,lplot,title5,x5,y5,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x6,y6,txres_6)
   dum7 = gsn_add_text(wks,lplot,title7,x7,y7,txres_7)
   dum8 = gsn_add_text(wks,lplot,title8,x8,y8,txres_8)
   dum9 = gsn_add_text(wks,lplot,title9,x9,y9,txres_9)
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

delete([/plot,dum0,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9/])
print("Figure done, yay!")

;;;;;;;;;;;;;
; Make difference plots
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_mass_budget_diffs"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

if (hemi .eq. "NH")
   pref = "d) Sea Ice Mass Budget Difference"
   ymin = -0.2
   ymax = 0.2
   yminor = 4
   x0 = 6.5 ;3.5
   x1 = 10.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 2.5 
   x6 = x5
   x7 = x5
   x8 = x1
   x9 = x5
   y0 = 0.18
   y1 = -0.06 ;0.1
   y2 = -0.08 ;0.08
   y3 = -0.1  ;0.06
   y4 = -0.12 ;0.04
   y5 = -0.06 ;0.12
   y6 = -0.08 ;0.1
   y7 = -0.1  ;0.08
   y8 = -0.14 ;0.04
   y9 = -0.12 ;0.06
end if
if (hemi .eq. "SH")
   pref = "b) Mass Budget Difference"
   ymin = -0.015
   ymax = 0.015
   yminor = 3
   x0 = 6.5
   x1 = 7.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 10.0
   x6 = x5
   x7 = x5
   x8 = x5
   x9 = x5
   y0 = 0.013
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

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,diff_plots_wrap,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref,x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   dum5 = gsn_add_text(wks,lplot,title5,x5,y5,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x6,y6,txres_6)
   dum7 = gsn_add_text(wks,lplot,title7,x7,y7,txres_7)
   dum8 = gsn_add_text(wks,lplot,title8,x8,y8,txres_8)
   dum9 = gsn_add_text(wks,lplot,title9,x9,y9,txres_9)
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

;;;;;;;;;;;;;
; Make difference plots - just biggest terms
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_mass_budget_diffs_subset"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

  pref = "b) Sea Ice Mass Budget Difference"
  delete(opts_b@xyLineColors)
  delete(opts_b@xyMarkerColors)
  delete(opts_b@xyMarkers)
  delete(opts_b@xyLineThicknesses)
  ; change colors
  opts_b@xyLineColors = (/"black","blue","deepskyblue","firebrick","palevioletred1","black","blue","deepskyblue","firebrick","palevioletred1"/)
  opts_b@xyMarkerColors = opts_b@xyLineColors
  opts_b@xyMarkers = (/4,4,4,4,4,16,16,16,16,16/)
  opts_b@xyLineThicknesses = (/5.,4.,4.,4.,4.,1.,1.,1.,1.,1./)

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,diff_plots_sub,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref,x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum5 = gsn_add_text(wks,lplot,title5,x4,y4,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x8,y8,txres_6)
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

;;;;;;;;;;;;;;;;;;;;;; END script
end
