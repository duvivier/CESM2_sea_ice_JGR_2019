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
; load albedo difference from atmosphere
; over: ice >= 15% concentration + open water
f1 = addfile(data_dir+"atmosphere_energy_budget_annual_cycle/CMIP6_b_bw_"+hemi+"_piControl_energy_budget_terms_atm_for_plot.nc","r")
b_data_all_atm = f1->b_data_all
bw_data_all_atm = f1->bw_data_all
diff_plot_atm = f1->diff_plots
delete(f1)

; load albedo difference from ice
; over: ice >= 15% concentration only
; (DOES NOT INCLUDE IMPACT FROM ICE FRACTION CHANGE)
f2 = addfile(data_dir+"ice_energy_budget_annual_cycle/CMIP6_b_bw_"+hemi+"_piControl_energy_budget_terms_for_plot.nc","r")
b_data_all_ice = f2->b_data_all
bw_data_all_ice = f2->bw_data_all
diff_plot_ice = f2->diff_plots
delete(f2)

; order in the files:
; 0=net; 1=swdn; 2=swup; 3=lwdn; 4=lwup; 5=conductive flux; 6=latent; 7=sensible

;;;;;;;;;;
; Combine for plotting
;;;;;;;;;;
nmons = 12
; For the plots we just want to plot the net and radiative (SW,LW) differences
diff_all = new((/20,nmons+1/),float)
; ice only will always be first
; net sfc energy
diff_all(0,0:nmons-1) = diff_plot_ice(0,:)  ; absolute
diff_all(1,0:nmons-1) = diff_plot_ice(8,:)  ; significance
diff_all(2,0:nmons-1) = diff_plot_atm(0,:)  ; absolute
diff_all(3,0:nmons-1) = diff_plot_atm(8,:)  ; significance
; swdn
diff_all(4,0:nmons-1) = diff_plot_ice(1,:)  ; absolute
diff_all(5,0:nmons-1) = diff_plot_ice(9,:)  ; significance
diff_all(6,0:nmons-1) = diff_plot_atm(1,:)  ; absolute
diff_all(7,0:nmons-1) = diff_plot_atm(9,:)  ; significance
; swup
diff_all(8,0:nmons-1) = diff_plot_ice(2,:)  ; absolute
diff_all(9,0:nmons-1) = diff_plot_ice(10,:)  ; significance
diff_all(10,0:nmons-1) = diff_plot_atm(2,:)  ; absolute
diff_all(11,0:nmons-1) = diff_plot_atm(10,:)  ; significance
; lwdn
diff_all(12,0:nmons-1) = diff_plot_ice(3,:)  ; absolute
diff_all(13,0:nmons-1) = diff_plot_ice(11,:)  ; significance
diff_all(14,0:nmons-1) = diff_plot_atm(3,:)  ; absolute
diff_all(15,0:nmons-1) = diff_plot_atm(11,:)  ; significance
; lwup
diff_all(16,0:nmons-1) = diff_plot_ice(4,:)  ; absolute
diff_all(17,0:nmons-1) = diff_plot_ice(12,:)  ; significance
diff_all(18,0:nmons-1) = diff_plot_atm(4,:)  ; absolute
diff_all(19,0:nmons-1) = diff_plot_atm(12,:)  ; significance

; now loop back to include January again
; ice only will always be first
; net sfc energy
diff_all(0,nmons) = diff_plot_ice(0,0)  ; absolute
diff_all(1,nmons) = diff_plot_ice(8,0)  ; significance
diff_all(2,nmons) = diff_plot_atm(0,0)  ; absolute
diff_all(3,nmons) = diff_plot_atm(8,0)  ; significance
; swdn
diff_all(4,nmons) = diff_plot_ice(1,0)  ; absolute
diff_all(5,nmons) = diff_plot_ice(9,0)  ; significance
diff_all(6,nmons) = diff_plot_atm(1,0)  ; absolute
diff_all(7,nmons) = diff_plot_atm(9,0)  ; significance
; swup
diff_all(8,nmons) = diff_plot_ice(2,0)  ; absolute
diff_all(9,nmons) = diff_plot_ice(10,0)  ; significance
diff_all(10,nmons) = diff_plot_atm(2,0)  ; absolute
diff_all(11,nmons) = diff_plot_atm(10,0)  ; significance
; lwdn
diff_all(12,nmons) = diff_plot_ice(3,0)  ; absolute
diff_all(13,nmons) = diff_plot_ice(11,0)  ; significance
diff_all(14,nmons) = diff_plot_atm(3,0)  ; absolute
diff_all(15,nmons) = diff_plot_atm(11,0)  ; significance
; lwup
diff_all(16,nmons) = diff_plot_ice(4,0)  ; absolute
diff_all(17,nmons) = diff_plot_ice(12,0)  ; significance
diff_all(18,nmons) = diff_plot_atm(4,0)  ; absolute
diff_all(19,nmons) = diff_plot_atm(12,0)  ; significance

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
if (hemi .eq. "NH")
   pref = "c) Surface Energy Budget Difference" 
   ymin = -10.0
   ymax = 25.0 
   yminor = 4
   x0 = 6.5
   x1 = 3.0
   x2 = x1
   x3 = 10.5
   x4 = x3
   x5 = x3
   x6 = x3
   x7 = x3
   x8 = x3
   x9 = x3
   x10 = x3
   y0 = 23.5 
   y1 = 18.75
   y2 = 16.75
   y3 = 20.75
   y4 = 18.75
   y5 = 16.75
   y6 = 14.75
   y7 = 12.75
   y8 = 10.75
   y9 = 8.75
   y10 = 6.75
end if
if (hemi .eq. "SH")
   pref = "d) Energy Budget Difference" ;"b) Antarctic Energy Budget"
   ymin = -3.0
   ymax = 5.0
   yminor = 1
   x0 = 6.0 ;3.5
   x1 = 7.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 10.0
   x6 = x5
   x7 = x5
   x8 = x5
   x9 = x5
   x10 = x5
   y0 = 4.0
   y1 = 3.75
   y2 = 3.25
   y3 = 2.75
   y4 = 2.25
   y5 = 3.75
   y6 = 3.25
   y7 = 2.75
   y8 = 2.25
   y9 = 2.25
   y10 = 2.25
end if

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
  opts@trYMinF                     = ymin
  opts@trYMaxF                     = ymax
  opts@tmYLMinorPerMajor           = yminor
  opts@tmYRMinorPerMajor           = yminor

; line resources for all ensembles
  opts@xyMarkerSizeF               = 0.025
  opts@tiMainOn = False
  opts@xyMarkLineMode = "MarkLines"
  opts@xyLineColors = (/"black","black","grey50","grey50","firebrick","firebrick","palevioletred1","palevioletred1","darkorange","darkorange","gold","gold","blue","blue","deepskyblue","deepskyblue","purple4","purple4","mediumpurple2","mediumpurple2"/)
  opts@xyMarkerColors = opts@xyLineColors
  opts@xyMarkers = (/4,16,4,16,4,16,4,16,4,16,4,16,4,16,4,16,4,16,4,16/)
  opts@xyLineThicknesses = (/4,0.1,4,0.1,4,0.1,4,0.1,4,0.1,4,0.1,4,0.1,4,0.1,4,0.1,4,0.1/)
  opts@xyDashPatterns = (/0,0,2,2,0,0,2,2,0,0,2,2,0,0,2,2,0,0,2,2/)

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2

; text
  txres_0 = txres
  txres_0@txFontHeightF = 0.03
  txres_0@txFontColor = "black"
  title1 = "Net (ice only)"
  txres_1 = txres
  txres_1@txFontColor = "black"
  title2 = "Net (ice + open water)"
  txres_2 = txres
  txres_2@txFontColor = "grey50"
  title3 = "SWdn (ice only)"
  txres_3 = txres
  txres_3@txFontColor = "firebrick"
  title4 = "SWdn (ice + open water)"
  txres_4 = txres
  txres_4@txFontColor = "palevioletred1"
  title5 = "SWup (ice only)"
  txres_5 = txres
  txres_5@txFontColor = "darkorange"
  title6 = "SWup (ice + open water)"
  txres_6 = txres
  txres_6@txFontColor = "gold"
  title7 = "LWdn (ice only)"
  txres_7 = txres
  txres_7@txFontColor = "blue"
  title8 = "LWdn (ice + open water)"
  txres_8 = txres
  txres_8@txFontColor = "deepskyblue"
  title9 = "LWup (ice only)"
  txres_9 = txres
  txres_9@txFontColor = "purple4"
  title10 = "LWup (ice + open water)"
  txres_10 = txres
  txres_10@txFontColor = "mediumpurple2"

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
; Make difference plots
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_energy_budget_diffs"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

   ; do actual plotting
   lplot = gsn_csm_xy(wks,xdata,diff_all,opts)
   dum0 = gsn_add_text(wks,lplot,pref,x0,y0,txres_0)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   dum5 = gsn_add_text(wks,lplot,title5,x5,y5,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x6,y6,txres_6)
   dum7 = gsn_add_text(wks,lplot,title7,x7,y7,txres_7)
   dum8 = gsn_add_text(wks,lplot,title8,x8,y8,txres_8)
   dum9 = gsn_add_text(wks,lplot,title9,x9,y9,txres_9)
   dum10 = gsn_add_text(wks,lplot,title10,x10,y10,txres_10)
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
