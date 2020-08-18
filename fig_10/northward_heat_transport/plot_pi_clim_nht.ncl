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
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
;***************************************
begin
;;;;;;;;;;
; Manual inputs or from cshell script
;;;;;;;;;;
; ******************************************************
; BEGIN SCRIPT

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
;data_dir = "./"
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_ice/analysis/version_9/northward_heat_transport/"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; load albedo difference from atmosphere
; over: ice >= 15% concentration + open water
f1 = addfile(data_dir+"CMIP6_b_bw_northward_heat_transport.nc","r")
b_data_all = f1->b_data_all
bw_data_all = f1->bw_data_all
diff = f1->diff
b_data_all_wm2 = f1->b_data_all_wm2
bw_data_all_wm2 = f1->bw_data_all_wm2
diff_wm2 = f1->diff_wm2
delete(f1)

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
cmip_dir = "/glade/collections/cdg/timeseries-cmip6/"
case1 = "b.e21.B1850.f09_g17.CMIP6-piControl.001"
b_files1 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*.FSNT.*.nc")
; get atmosphere dimensions
fin = addfile(b_files1(0),"r")
lat1d_atm = fin->lat
lon1d_atm = fin->lon
gw = tofloat(fin->gw)
nlat_atm = dimsizes(lat1d_atm)
nlon_atm = dimsizes(lon1d_atm)
delete(fin)

;;;;;;;;;;
; Set constants
;;;;;;;;;;
    pi      = 3.14159265    ; pi
    g       = 9.81          ; gravity (m/s2)
    Lf      = 3.337e5       ; latent heat of fusion (J/kg)
    Re      = 6.371e6       ; radius of earth (m)
    rho_h20 = 1000.         ; density of water (kg/m3)
    rho_ice = 917.          ; density of ice (kg/m3)

    ; derived values for NHT
    dlon = (2.*pi)/nlon_atm          ; change in longitude (radians) per unit latitude
    coef = (Re^2)/1.e15             ; (m2 and scaled for PW)
    ; derived values for ice NHT
    Earea = 4.*pi*Re^2               ; surface area of Earth (m2)
    dlon_deg = (lon1d_atm(2) - lon1d_atm(1))  ; assume dlon is constant
    dlat_deg = (lat1d_atm(2) - lon1d_atm(1))  ; assume dlat is constant
    wgt = tofloat(NormCosWgtGlobe(lat1d_atm))
    tsum = sum(wgt)          ; sum of all weights
    nwgt = wgt/tsum          ; fraction of sphre of each latitude band
    boxarea = Earea*nwgt/nlon_atm  ; area of each grid box latitude in m2

    ; get surface area of Earth above each latitude band
    ; A = 2*pi*Re^2*(1-cos(90-lat))  - NOTE, must convert latitude to radians
    lat_in = tofloat((90.0-lat1d_atm)*(pi/180.))  ; get latitude to input in radians
    area_sfc_cap = (2*pi*(Re^2))*(1-cos(lat_in))  ; in m2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting NHT")
; set outdir and fout names
   outdir = "./"

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
  opts@tiMainOn                    = False
  opts@tmXBOn                      = True
  opts@trXMinF                     = -90.0
  opts@trXMaxF                     = 90.0
  opts@tiXAxisString               = "Latitude"
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@gsnYRefLine                 = 0.0

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025
  opts@xyMarkLineMode              = "Lines"
  opts@xyLineColors                = (/"black","firebrick","blue","deepskyblue"/)
  opts@xyLineThicknesses           = (/5.0,4.0,4.0,4.0/)
  opts@xyDashPatterns              = (/0,0,0,0/)

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2

; text
  title1 = "Net NHT"
  txres_1 = txres
  txres_1@txFontColor = "black"
  title2 = "ATM NHT"
  txres_2 = txres
  txres_2@txFontColor = "firebrick"
  title3 = "OCN NHT"
  txres_3 = txres
  txres_3@txFontColor = "blue"
  title4 = "ICE NHT"
  txres_4 = txres
  txres_4@txFontColor = "deepskyblue"

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 3                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 3                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True	       ; fill the page

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Make plots for NHT in PW
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
opts_b = opts
units = "PW"
;;;;;;;;;;;;;
; Make b case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_full_piControl_NHT"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)  

   ymin = -6.0
   ymax = 6.0
   yminor = 3
   x0 = -65.0
   x1 = -60.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 5.5
   y1 = 4.0
   y2 = 3.25
   y3 = 2.5
   y4 = 1.75 

  opts_b@tiYAxisString               = units
  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor  

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,b_data_all,opts_b)
   dum0 = gsn_add_text(wks,lplot,"a) NHT CAM6",x0,y0,txres_1)
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

;;;;;;;;;;;;;
;  Make bw case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_bw_full_piControl_NHT"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)    

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,bw_data_all,opts_b)
   dum0 = gsn_add_text(wks,lplot,"a) NHT WACCM6",x0,y0,txres_1)
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

;;;;;;;;;;;;;
; Make difference plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_full_piControl_NHT_diffs"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

   ymin = -0.1
   ymax = 0.1
   yminor = 3
   x0 = 0.0 ;-60.0
   x1 = -40.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 0.085
   y1 = -0.04
   y2 = -0.05
   y3 = -0.06
   y4 = -0.07

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor   

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,diff,opts_b)
   dum0 = gsn_add_text(wks,lplot,"b) NHT Difference",x0,y0,txres_1)
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

;;;;;;;;;;;;;
; Make difference plot in Arctic only
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_full_piControl_NHT_Arctic_diffs"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

   ymin = -0.08
   ymax = 0.02
   yminor = 3
   x0 = 70.0
   x1 = 80.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 0.015
   y1 = -0.04
   y2 = -0.05
   y3 = -0.06
   y4 = -0.07

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor
  opts_b@trXMinF                     = 55.0
  opts_b@trXMaxF                     = 85.0

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,diff,opts_b)
   dum0 = gsn_add_text(wks,lplot,"b) NHT Difference",x0,y0,txres_1)
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Make plots for NHT in W/m2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
delete(opts_b)
opts_b = opts
units = "W m~S~-2~N~"
;;;;;;;;;;;;;
; Make b case plot 
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_full_piControl_NHT_wm2"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)  

   ymin = -50.0
   ymax = 300.0
   yminor = 4
   x0 = -65.0
   x1 = -60.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 275.0
   y1 = 175.0
   y2 = 150.0
   y3 = 125.0
   y4 = 100.0

   opts_b@tiYAxisString               = units
   opts_b@trYMinF                     = ymin
   opts_b@trYMaxF                     = ymax
   opts_b@tmYLMinorPerMajor           = yminor
   opts_b@tmYRMinorPerMajor           = yminor  

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,b_data_all_wm2,opts_b)
   dum0 = gsn_add_text(wks,lplot,"a) NHT CAM6",x0,y0,txres_1)
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

;;;;;;;;;;;;;
;  Make bw case plot - NHT in PW
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_bw_full_piControl_NHT_wm2"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)    

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,bw_data_all_wm2,opts_b)
   dum0 = gsn_add_text(wks,lplot,"a) NHT WACCM6",x0,y0,txres_1)
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

;;;;;;;;;;;;;
; Make difference plot - in PW
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_full_piControl_NHT_diffs_wm2"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

   ymin = -8.0
   ymax = 2.0
   yminor = 3
   x0 = 0.0 
   x1 = -40.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 1.5
   y1 = -2.0
   y2 = -3.0
   y3 = -4.0
   y4 = -5.0

   opts_b@trYMinF                     = ymin
   opts_b@trYMaxF                     = ymax
   opts_b@tmYLMinorPerMajor           = yminor
   opts_b@tmYRMinorPerMajor           = yminor   

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,diff_wm2,opts_b)
   dum0 = gsn_add_text(wks,lplot,"b) NHT Difference",x0,y0,txres_1)
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

;;;;;;;;;;;;;
; Make difference plot in Arctic only - in PW
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_full_piControl_NHT_Arctic_diffs_wm2"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

   ymin = -5.0
   ymax = 1.5
   yminor = 3
   x0 = 70.0
   x1 = 60.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 1.0
   y1 = -2.5
   y2 = -3.0
   y3 = -3.5
   y4 = -4.0

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor
  opts_b@trXMinF                     = 55.0
  opts_b@trXMaxF                     = 85.0   

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,diff_wm2,opts_b)
   dum0 = gsn_add_text(wks,lplot,"a) NHT Difference",x0,y0,txres_1)
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

;;;;;;;;;;;;;;;;;;;;;; END script
end
