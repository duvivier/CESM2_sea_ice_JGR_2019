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
varcode = "sfc_snow"
hemi = "NH"     ; choose hemisphere - "NH" or "SH"
; ******************************************************
; BEGIN SCRIPT

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_ice/analysis/version_9/snow_annual_cycle/"

nmons = 12
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; snow mass budget terms
f1 = addfile(data_dir+"/snow_mass_budget_annual_cycle/CMIP6_b_bw_"+hemi+"_piControl_snow_mass_budget_terms_for_plot.nc","r")
diff_plots = f1->diff_plots
diff_melt_rate = new((/2,nmons+1/),float)
diff_melt_rate(0,0:nmons-1) = diff_plots(2,:)
diff_melt_rate(1,0:nmons-1) = diff_plots(6,:)
diff_melt_rate(0,nmons) = diff_plots(2,0)
diff_melt_rate(1,nmons) = diff_plots(6,0)
delete(f1)
; snow depth
f2 = addfile(data_dir+"/snow_depth_annual_cycle/CMIP6_b_bw_"+hemi+"_piControl_sisnthick_diff_annual_cycle_for_plot.nc","r")
diff = f2->diff
; convert m to cm and wrap months
diff_depth = new((/2,nmons+1/),float)
diff_depth(:,0:nmons-1) = diff*100.
diff_depth(:,nmons) = diff(:,0)*100.
delete(f2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"

xdata = ispan(1, nmons+1, 1)
label_inds = xdata
label_str = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)

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
  opts@tmXTOn                      = True
  opts@tmXTMode                    = "Explicit"
  opts@tmXTValues                  = label_inds
  opts@tmXTLabels                  = label_str
  opts@trXMinF                     = 1
  opts@trXMaxF                     = nmons+1
  opts@tiXAxisString               = "Month"
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tmXTLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@tmYRLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@gsnYRefLine                 = 0.0

; resources for budget terms
  opts_b = opts
  opts_b@tiYAxisString               = "cm"
  opts_b@trYMinF                     = -2.0
  opts_b@trYMaxF                     = 2.0
  opts_b@tiMainOn                    = False
  opts_b@xyMarkLineMode              = "MarkLines"
  opts_b@xyLineColors                = (/"blue","blue"/)
  opts_b@xyMarkerColors              = opts_b@xyLineColors
  opts_b@xyMarkers                   = (/4,16/)
  opts_b@xyMarkerSizeF               = 0.025
  opts_b@xyLineThicknesses           = (/4.,4./)
  opts_b@xyDashPattern               = 0
; resources for snow_depth
  opts_bw = opts
  opts_bw@tiYAxisString               = "cm/day"
  opts_bw@trYMinF                     = -0.035
  opts_bw@trYMaxF                     = 0.035
  opts_bw@tiMainOn                    = False
  opts_bw@xyMarkLineMode              = "MarkLines"
  opts_bw@xyLineColors                = (/"firebrick","firebrick"/)
  opts_bw@xyMarkerColors              = opts_bw@xyLineColors
  opts_bw@xyMarkers                   = (/4,16/)
  opts_bw@xyMarkerSizeF               = 0.025
  opts_bw@xyLineThicknesses           = (/4.,4./)
  opts_bw@xyDashPattern               = 0

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2
  title1 = "c) Snow on Sea Ice Difference"
  txres_1 = txres
  txres_1@txFontColor = "black"
  x1 = 6.0
  y1 = 0.03
  title2 = "Snow Melt Rate (cm/day)"
  txres_2 = txres
  x2 = 3.5
  y2 = 0.02
  txres_2@txFontColor = "firebrick"
  title3 = "Snow Depth (cm)"
  txres_3 = txres
  txres_3@txFontColor = "blue"
  x3 = 3.5
  y3 = 0.015

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
; Make b case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_"+hemi+"_piControl_"+varcode
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic) 

   ; do actual plotting
   lplot = gsn_csm_xy2(wks,xdata,diff_melt_rate,diff_depth,opts_bw,opts_b)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
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
