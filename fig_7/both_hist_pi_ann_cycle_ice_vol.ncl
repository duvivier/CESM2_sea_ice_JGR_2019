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
; Manual inputs
;;;;;;;;;;
varcode = "sivol"
hemi = "SH"     ; choose hemisphere - "NH" or "SH"

;;;;;;;;;;
; Set some info for plots etc.
;;;;;;;;;;
varcode2 = "siconc"
vartitle = "ice_vol"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "./"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; load historical
f1 = addfile(data_dir+"CMIP6_b_bw_"+hemi+"_hist_"+vartitle+"_for_annual_cycle.nc","r")
; load data from historical file
all_vol_hemi_mon = f1->all_vol_hemi_mon
all_vol_hemi_stddev = f1->all_vol_hemi_stddev
delete(f1)
; 0=b; 1=bw; 2=nsidc, 3=masked_b, 4=masked_bw

; load pi
f2 = addfile(data_dir+"CMIP6_b_bw_"+hemi+"_piControl_"+vartitle+"_for_annual_cycle.nc","r")
; load data from historical file
data_sub = f2->data_sub
data_stddev = f2->data_stddev
delete(f2)
; 0=b; 1=bw; 2=masked_b; 3=masked_bw

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Combine to one array
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nmons = 12
mean_all = new((/10,nmons+1/),float)
mean_all(0:3,0:nmons-1) = data_sub
mean_all(0:3,nmons) = data_sub(:,0)
mean_all(4:5,0:nmons-1) = all_vol_hemi_mon(0:1,:)
mean_all(4:5,nmons) = all_vol_hemi_mon(0:1,0)
mean_all(6:7,0:nmons-1) = all_vol_hemi_mon(4:5,:)
mean_all(6:7,nmons) = all_vol_hemi_mon(4:5,0)
mean_all(8:9,0:nmons-1) = all_vol_hemi_mon(2:3,:)
mean_all(8:9,nmons) = all_vol_hemi_mon(2:3,0)

stddev_all = new((/10,nmons+1/),float)
stddev_all(0:3,0:nmons-1) = data_stddev
stddev_all(0:3,nmons) = data_stddev(:,0)
stddev_all(4:5,0:nmons-1) = all_vol_hemi_stddev(0:1,:)
stddev_all(4:5,nmons) = all_vol_hemi_stddev(0:1,0)
stddev_all(6:7,0:nmons-1) = all_vol_hemi_stddev(4:5,:)
stddev_all(6:7,nmons) = all_vol_hemi_stddev(4:5,0)
stddev_all(8:9,0:nmons-1) = all_vol_hemi_stddev(2:3,:)
stddev_all(8:9,nmons) = all_vol_hemi_stddev(2:3,0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_"+vartitle+"_annual_cycle"

xdata = ispan(1, 13, 1)
label_inds = xdata
label_str = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)

if (hemi .eq. "NH")
   pref0 = "c) " ;Arctic - "
   pref1 = "d) " ;Arctic - "
   ymin = 5.0
   ymax = 38.0
   yminor = 4
   ymin1 = 1.5
   ymax1 = 6.0
   yminor1 = 3
   x0 = 6.0
   x1 = 4.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = x1
   x6 = x1
   x7 = 6.0
   y0 = 36.5
   y1 = 17.0
   y2 = 15.0
   y3 = 13.0
   y4 = 11.0
   y5 = 9.0
   y6 = 7.0
   y7 = 5.75
end if
if (hemi .eq. "SH")
   pref0 = "c) " ;Antarctic - "
   pref1 = "d) " ;Antarctic - "
   ymin = 0.0
   ymax = 24.0
   yminor = 3
   ymin1 = 0.4
   ymax1 = 1.4
   yminor1 = 3
   x0 = 6.0 
   x1 = 9.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = x1
   x6 = x1
   x7 = 6.0
   y0 = 22.5
   y1 = 12.0
   y2 = 10.5
   y3 = 9.0
   y4 = 7.5
   y5 = 6.0
   y6 = 4.5
   y7 = 1.35
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
  opts@tiYAxisString               = "10~S~3~N~ km~S~3~N~"
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@trYMinF                     = ymin
  opts@trYMaxF                     = ymax
  opts@tmYLMinorPerMajor           = yminor
  opts@tmYRMinorPerMajor           = yminor

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2

; resources for nsidc
  opts = opts
  opts@tiMainOn = False
  opts@xyMarkLineMode = "MarkLines"
  ; black = CAM historical; firebrick = WACCM historical; grey50 = CAM PI; lightpink3 = WACCM PI, blue = NSIDC
  opts@xyLineColors = (/"grey50","lightpink3","grey50","lightpink3","black","firebrick","black","firebrick","blue","deepskyblue"/)
  opts@xyLineThicknesses = (/3,3,0.01,0.01,3,3,0.01,0.01,3,3/)
  opts@xyDashPatterns = (/2,2,2,2,0,0,0,0,0,0/)
  opts@xyMarkers = (/4,4,16,16,4,4,16,16,1,1/) 
  opts@xyMarkerColors = opts@xyLineColors

; text resources
  txres_o = txres
  txres_o@txFontColor = "blue"
  txres_o2 = txres
  txres_o2@txFontColor = "deepskyblue"
  txres_bp = txres
  txres_bp@txFontColor = "grey50"
  txres_bh = txres
  txres_bh@txFontColor = "black"
  txres_bwp = txres
  txres_bwp@txFontColor = "lightpink3"
  txres_bwh = txres
  txres_bwh@txFontColor = "firebrick"

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 3                ; Add white space b/w plos.
  pnlres@gsnPanelXWhiteSpacePercent = 3                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True	       ; fill the page

;;;;;;;;;;;;;
; Make plots
;;;;;;;;;;;;;
 wks = gsn_open_wks("ps" ,fout)
 plot = new(2, graphic)    

   ; do actual plotting for both seasonal cycles
   p = 0
   lplot = gsn_csm_xy(wks,xdata,mean_all,opts)
   ;dum0 = gsn_add_text(wks,lplot,pref0+"Mean Ice Volume",x0,y0,txres_bh)
   dum1 = gsn_add_text(wks,lplot,"GIOMAS",x1,y1,txres_o)
   dum2 = gsn_add_text(wks,lplot,"PIOMAS",x2,y2,txres_o2)
   dum3 = gsn_add_text(wks,lplot,"CAM6 - historical",x3,y3,txres_bh)
   dum4 = gsn_add_text(wks,lplot,"WACCM6 - historical",x4,y4,txres_bwh)
   dum5 = gsn_add_text(wks,lplot,"CAM6 - PI",x5,y5,txres_bp)
   dum6 = gsn_add_text(wks,lplot,"WACCM6 - PI",x6,y6,txres_bwp)
   plot(p) = lplot

   ;seasonal stddev
   opts@trYMinF                     = ymin1
   opts@trYMaxF                     = ymax1
   opts@tmYLMinorPerMajor           = yminor1
   opts@tmYRMinorPerMajor           = yminor1
   p = 1
   lplot = gsn_csm_xy(wks,xdata,stddev_all,opts)
   ;dum7 = gsn_add_text(wks,lplot,pref1+"StdDev. Ice Volume",x7,y7,txres_bh)
   plot(p) = lplot

;;;;;;;;;;;;;
; Final Plots
;;;;;;;;;;;;;
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,2/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Plot 1 complete")

;;;;;;;;;;;;;;;;;;;;;; END script
end
