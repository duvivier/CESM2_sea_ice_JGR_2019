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
varcode = "siconc"
hemi = "NH"     ; choose hemisphere - "NH" or "SH"

;;;;;;;;;;
; Set some info for plots etc.
;;;;;;;;;;
varcode2 = "siconc"
vartitle = "ice_ext"

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
all_ext_hemi_mon = f1->all_ext_hemi_mon
all_ext_hemi_stddev = f1->all_ext_hemi_stddev
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
; Combine to one array and wrap around back to Jan
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nmons = 12
mean_all = new((/9,nmons+1/),float)
mean_all(0:3,0:nmons-1) = data_sub
mean_all(0:3,nmons) = data_sub(:,0)
mean_all(4:5,0:nmons-1) = all_ext_hemi_mon(0:1,:)
mean_all(4:5,nmons) = all_ext_hemi_mon(0:1,0)
mean_all(6:7,0:nmons-1) = all_ext_hemi_mon(3:4,:)
mean_all(6:7,nmons) = all_ext_hemi_mon(3:4,0)
mean_all(8,0:nmons-1) = all_ext_hemi_mon(2,:)
mean_all(8,nmons) = all_ext_hemi_mon(2,0)

stddev_all = new((/9,nmons+1/),float)
stddev_all(0:3,0:nmons-1) = data_stddev
stddev_all(0:3,nmons) = data_stddev(:,0)
stddev_all(4:5,0:nmons-1) = all_ext_hemi_stddev(0:1,:)
stddev_all(4:5,nmons) = all_ext_hemi_stddev(0:1,0)
stddev_all(6:7,0:nmons-1) = all_ext_hemi_stddev(3:4,:)
stddev_all(6:7,nmons) = all_ext_hemi_stddev(3:4,0)
stddev_all(8,0:nmons-1) = all_ext_hemi_stddev(2,:)
stddev_all(8,nmons) = all_ext_hemi_stddev(2,0)

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
   pref0 = "a) " ;Arctic - "
   pref1 = "b) " ;Arctic - "
   ymin = 4.0
   ymax = 18.0
   yminor = 1
   ymin1 = 0.0
   ymax1 = 1.8
   yminor1 = 2
   x0 = 6.0
   x1 = 3.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = x1
   x6 = 6.0
   y0 = 17.0
   y1 = 12.0
   y2 = 11.0
   y3 = 10.0
   y4 = 9.0
   y5 = 8.0
   y6 = 1.65
end if
if (hemi .eq. "SH")
   pref0 = "a) " ;Antarctic - "
   pref1 = "b) " ;Antarctic - "
   ymin = 0.0
   ymax = 20.0
   yminor = 2
   ymin1 = 0.3
   ymax1 = 1.0
   yminor1 = 1
   x0 = 6.0
   x1 = 9.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = x1
   x6 = 6.0
   y0 = 19.0
   y1 = 8.0
   y2 = 7.0
   y3 = 6.0
   y4 = 5.0
   y5 = 4.0
   y6 = 0.95
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
  opts@tiYAxisString               = "10~S~6~N~ km~S~2~N~"
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
  opts@xyLineColors = (/"grey50","lightpink3","grey50","lightpink3","black","firebrick","black","firebrick","blue"/)
  opts@xyLineThicknesses = (/3,3,0.01,0.01,3,3,0.01,0.01,3/)
  opts@xyDashPatterns = (/2,2,2,2,0,0,0,0,0/)
  opts@xyMarkers = (/4,4,16,16,4,4,16,16,1/) 
  opts@xyMarkerColors = opts@xyLineColors

; text resources
  txres_o = txres
  txres_o@txFontColor = "blue"
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
   ;dum0 = gsn_add_text(wks,lplot,pref0+"Mean Ice Extent",x0,y0,txres_bh)
   dum1 = gsn_add_text(wks,lplot,"NSIDC",x1,y1,txres_o)
   dum2 = gsn_add_text(wks,lplot,"CAM6 - historical",x2,y2,txres_bh)
   dum3 = gsn_add_text(wks,lplot,"WACCM6 - historical",x3,y3,txres_bwh)
   dum4 = gsn_add_text(wks,lplot,"CAM6 - PI",x4,y4,txres_bp)
   dum5 = gsn_add_text(wks,lplot,"WACCM6 - PI",x5,y5,txres_bwp)
   plot(p) = lplot

   ;seasonal stddev
   opts@trYMinF                     = ymin1
   opts@trYMaxF                     = ymax1
   opts@tmYLMinorPerMajor           = yminor1
   opts@tmYRMinorPerMajor           = yminor1
   p = 1
   lplot = gsn_csm_xy(wks,xdata,stddev_all,opts)
   ;dum6 = gsn_add_text(wks,lplot,pref1+"StdDev. Ice Extent",x6,y6,txres_bh)
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
