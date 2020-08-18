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
varcode = "net_energy_budget"
hemi = "SH"     ; choose hemisphere - "NH" or "SH"
; ******************************************************
; BEGIN SCRIPT

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
;data_dir = "./"
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_ice/analysis/version_9/energy_budget_annual_cycle/"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; load albedo difference from ice
; over: ice >= 15% concentration only
; (DOES NOT INCLUDE IMPACT FROM ICE FRACTION CHANGE)
f2 = addfile(data_dir+"ice_energy_budget_annual_cycle/CMIP6_b_bw_"+hemi+"_piControl_energy_budget_terms_for_plot.nc","r")
b_data_all_in = f2->b_data_all
bw_data_all_in = f2->bw_data_all
diff_plots_in = f2->diff_plots
delete(f2)

; order in the files:
; 0=net; 1=swdn; 2=swup; 3=lwdn; 4=lwup; 5=conductive flux; 6=latent; 7=sensible

; now loop back to include January again
nmons = 12
b_data_all = new((/8,nmons+1/),float)
b_data_all(:,0:nmons-1) = b_data_all_in
b_data_all(:,nmons) = b_data_all_in(:,0)

bw_data_all = new((/8,nmons+1/),float)
bw_data_all(:,0:nmons-1) = bw_data_all_in
bw_data_all(:,nmons) = bw_data_all_in(:,0)

diff_plots = new((/16,nmons+1/),float)
diff_plots(:,0:nmons-1) = diff_plots_in
diff_plots(:,nmons) = diff_plots_in(:,0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"

xdata = ispan(1, 13, 1)
label_inds = xdata
label_str = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)

units = "W m~S~-2~N~"

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
  opts_b@xyLineColors = (/"black","firebrick","palevioletred1","darkorange","gold","green3","blue","deepskyblue","black","firebrick","palevioletred1","darkorange","gold","green3","blue","deepskyblue"/)
  opts_b@xyMarkerColors = opts_b@xyLineColors
  opts_b@xyMarkers = (/1,1,1,1,1,1,1,1,16,16,16,16,16,16,16,16/)
  opts_b@xyLineThicknesses = (/5.,4.,4.,4.,4.,4.,4.,4.,1.,1.,1.,1.,1.,1.,1.,1./)
  opts_b@xyDashPattern = 0

; text
  title1 = "Net"
  txres_1 = txres
  txres_1@txFontColor = "black"
  title2 = "Shortwave Down"
  txres_2 = txres
  txres_2@txFontColor = "firebrick"
  title3 = "Shortwave Up"
  txres_3 = txres
  txres_3@txFontColor = "palevioletred1"
  title4 = "Longwave Down"
  txres_4 = txres
  txres_4@txFontColor = "darkorange"
  title5 = "Longwave Up"
  txres_5 = txres
  txres_5@txFontColor = "gold"
  title6 = "Top Conduction"
  txres_6 = txres
  txres_6@txFontColor = "green3"
  title7 = "Latent Heat"
  txres_7 = txres
  txres_7@txFontColor = "blue"
  title8 = "Sensible Heat"
  txres_8 = txres
  txres_8@txFontColor = "deepskyblue"

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
   pref = "e) ";"b) Arctic "
   ymin = -100.0
   ymax = 500.0 ;400.0
   yminor = 3
   x0 = 3.5
   x1 = 8.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 11.
   x6 = x5
   x7 = x5
   x8 = x5
   y0 = 450.
   y1 = 475.
   y2 = 435.
   y3 = 395.
   y4 = 355.
   y5 = 475.
   y6 = 435.
   y7 = 395.
   y8 = 355.
end if
if (hemi .eq. "SH")
   pref = "c) " ; "a) Antarctic"
   ymin = -100.0
   ymax = 500.0
   yminor = 3
   x0 = 3.5
   x1 = 8.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 11.
   x6 = x5
   x7 = x5
   x8 = x5
   y0 = 450.
   y1 = 475.
   y2 = 435.
   y3 = 395.
   y4 = 355.
   y5 = 475.
   y6 = 435.
   y7 = 395.
   y8 = 355.
end if

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor

;;;;;;;;;;;;;
; Make b case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_"+hemi+"_piControl_energy_budget"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic) 

   ; reverse the swup and lwup for better axes   
   b_data_all(2,:) = -1.0*b_data_all(2,:)
   b_data_all(4,:) = -1.0*b_data_all(4,:)

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,b_data_all,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+" Sea Ice Energy Budget~C~CAM6",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,"-1*"+title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   dum5 = gsn_add_text(wks,lplot,"-1*"+title5,x5,y5,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x6,y6,txres_6)
   dum7 = gsn_add_text(wks,lplot,title7,x7,y7,txres_7)
   dum8 = gsn_add_text(wks,lplot,title8,x8,y8,txres_8)
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

delete([/plot,dum0,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8/])
print("Figure done, yay!")

;;;;;;;;;;;;;
; Make bw case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_bw_"+hemi+"_piControl_energy_budget"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic) 

   ; reverse the swup and lwup for better axes   
   bw_data_all(2,:) = -1.0*bw_data_all(2,:)
   bw_data_all(4,:) = -1.0*bw_data_all(4,:)

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,bw_data_all,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+" Sea Ice Energy Budget~C~WACCM6",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,"-1*"+title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   dum5 = gsn_add_text(wks,lplot,"-1*"+title5,x5,y5,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x6,y6,txres_6)
   dum7 = gsn_add_text(wks,lplot,title7,x7,y7,txres_7)
   dum8 = gsn_add_text(wks,lplot,title8,x8,y8,txres_8)
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

delete([/plot,dum0,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8/])
print("Figure done, yay!")

;;;;;;;;;;;;;
; Make difference plots
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_energy_budget_diffs"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

if (hemi .eq. "NH")
   pref = "f) Sea Ice Energy Budget Difference" ;"c) Arctic Energy Budget"
   ymin = -10.0
   ymax = 15.0
   yminor = 3
   x0 = 4.5
   x1 = 10.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = x1
   x6 = x1
   x7 = x1
   x8 = x1
   y0 = 13.0
   y1 = 14.0
   y2 = 12.5
   y3 = 11.
   y4 = 9.5
   y5 = 8.
   y6 = 6.5
   y7 = 5.
   y8 = 3.5
end if
if (hemi .eq. "SH")
   pref = "d) Energy Budget Difference" ;"b) Antarctic Energy Budget"
   ymin = -3.0
   ymax = 5.0
   yminor = 1
   x0 = 4.0
   x1 = 7.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 10.0
   x6 = x5
   x7 = x5
   x8 = x5
   y0 = 4.5
   y1 = 3.75
   y2 = 3.25
   y3 = 2.75
   y4 = 2.25
   y5 = 3.75
   y6 = 3.25
   y7 = 2.75
   y8 = 2.25
end if

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,diff_plots,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+"~C~(CAM6-WACCM6)",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   dum5 = gsn_add_text(wks,lplot,title5,x5,y5,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x6,y6,txres_6)
   dum7 = gsn_add_text(wks,lplot,title7,x7,y7,txres_7)
   dum8 = gsn_add_text(wks,lplot,title8,x8,y8,txres_8)
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
