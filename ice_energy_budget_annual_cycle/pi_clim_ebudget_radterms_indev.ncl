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
varcode = "silwnet"
hemi = "NH"     ; choose hemisphere - "NH" or "SH"
; ******************************************************
; BEGIN SCRIPT

if (varcode .eq. "siswnet")
   varcode1 = "siflswdtop"
   varcode2 = "siflswutop"
end if
if (varcode .eq. "silwnet")
   varcode1 = "sifllwdtop"
   varcode2 = "sifllwutop"
end if
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/NSIDC_SeaIce_extent/"
grid_dir = "/glade/p/cesm/omwg/grids/"

b_files1 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode1+"/gn/latest/*.nc")
b_files2 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode2+"/gn/latest/*.nc")
nfiles_b = dimsizes(b_files1)
bw_files1 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode1+"/gn/latest/*.nc")
bw_files2 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode2+"/gn/latest/*.nc")
nfiles_bw = dimsizes(bw_files1)

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(grid_dir+"gx1v7_grid.nc","r")
lat2d = b_file->TLAT
lon2d = b_file->TLONG
tarea = b_file->TAREA
tarea = tarea/100.0/100.0 ; convert from cm2 to m2
dims = dimsizes(lat2d)
nlat = dims(0)
nlon = dims(1)
delete(b_file)

; get hemispheric mask
mask_hemi = tofloat(lat2d)

if (hemi .eq. "NH")then
   mask_hemi = where(lat2d .ge. 30.0, 1.0, mask_hemi@_FillValue)
end if
if (hemi .eq. "SH")then
   mask_hemi = where(lat2d .le. -30.0, 1.0, mask_hemi@_FillValue)
end if

nmons = 12
nyrs_b = 1200
nyrs_bw = 499
ntimes_b = nmons*nyrs_b
ntimes_bw = nmons*nyrs_bw

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data, Calculate hemispheric means at the same time")
; data is in 50 or 100 year chunks
; avg over just hemi cells

;;; b files
print("Loading b expt PI control")
b_data = new((/ntimes_b/),"float")
b_data_1 = new((/ntimes_b/),"float")
b_data_2 = new((/ntimes_b/),"float")
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   fin1 = addfile(b_files1(n),"r")
   ; deal with variable of interest
   dum_1 = fin1->$varcode1$
   fin2 = addfile(b_files2(n),"r")
   ; deal with variable of interest
   dum_2 = fin2->$varcode2$
   if (varcode .eq. "siswnet")
      dum_2 = -1.0*dum_2
      dum_v = dum_1 + dum_2
   end if
   if (varcode .eq. "silwnet")
      dum_v = dum_1 + dum_2
   end if
   dum_1_hemi = dim_avg_n(dim_avg_n(dum_1*conform(dum_1,mask_hemi,(/1,2/)),2),1)
   dum_2_hemi = dim_avg_n(dim_avg_n(dum_2*conform(dum_2,mask_hemi,(/1,2/)),2),1)
   dum_v_hemi = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_hemi,(/1,2/)),2),1) 
   ; deal with time index
   dum_t = fin1->time
   ; deal with indices
   dim_t = dimsizes(dum_t)
   if (n .eq. 0)
      ind_st = n
      ind_ed = dim_t -1
   else
      ind_st = ind_ed + 1
      ind_ed = ind_st+dim_t - 1
   end if
   ; put data in the arrays
   b_data_1(ind_st:ind_ed) = dum_1_hemi
   b_data_2(ind_st:ind_ed) = dum_2_hemi
   b_data(ind_st:ind_ed) = dum_v_hemi
   b_time(ind_st:ind_ed) = dum_t
   delete([/fin1,fin2,dum_t,dum_1,dum_1_hemi,dum_2,dum_2_hemi,dum_v,dum_v_hemi,dim_t/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data = new((/ntimes_bw/),"float")
bw_data_1 = new((/ntimes_bw/),"float")
bw_data_2 = new((/ntimes_bw/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   fin1 = addfile(bw_files1(n),"r")
   ; deal with variable of interest
   dum_1 = fin1->$varcode1$
   fin2 = addfile(bw_files2(n),"r")
   ; deal with variable of interest
   dum_2 = fin2->$varcode2$
   if (varcode .eq. "siswnet")
      dum_2 = -1.0*dum_2
      dum_v = dum_1 + dum_2
   end if
   if (varcode .eq. "silwnet")
      dum_v = dum_1 + dum_2
   end if
   dum_1_hemi = dim_avg_n(dim_avg_n(dum_1*conform(dum_1,mask_hemi,(/1,2/)),2),1)
   dum_2_hemi = dim_avg_n(dim_avg_n(dum_2*conform(dum_2,mask_hemi,(/1,2/)),2),1)
   dum_v_hemi = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_hemi,(/1,2/)),2),1) 
   ; deal with time index
   dum_t = fin1->time
   ; deal with indices
   dim_t = dimsizes(dum_t)
   if (n .eq. 0)
      ind_st = n
      ind_ed = dim_t -1
   else
      ind_st = ind_ed + 1
      ind_ed = ind_st+dim_t - 1
   end if
   ; put data in the arrays
   b_data_1(ind_st:ind_ed) = dum_1_hemi
   b_data_2(ind_st:ind_ed) = dum_2_hemi
   bw_data(ind_st:ind_ed) = dum_v_hemi
   bw_time(ind_st:ind_ed) = dum_t
   delete([/fin1,fin2,dum_t,dum_1,dum_1_hemi,dum_2,dum_2_hemi,dum_v,dum_v_hemi,dim_t/])
end do

;;;;;;;;;;
; Convert time arrays
;;;;;;;;;;
b_time_conv = cd_calendar(b_time,0)
yy_b = b_time_conv(:,0)
mm_b = b_time_conv(:,1)

bw_time_conv = cd_calendar(bw_time,0)
yy_bw = bw_time_conv(:,0)
mm_bw = bw_time_conv(:,1)

;;;;;;;;;;
; Get climatologies
;;;;;;;;;;
print ("Getting annual climatology")
b_data_avg = new((/nmons/),typeof(b_data))
bw_data_avg = new((/nmons/),typeof(bw_data))
b_data_1_avg = new((/nmons/),typeof(b_data))
bw_data_1_avg = new((/nmons/),typeof(bw_data))
b_data_2_avg = new((/nmons/),typeof(b_data))
bw_data_2_avg = new((/nmons/),typeof(bw_data))
b_data_stddev = new((/nmons/),typeof(b_data))
bw_data_stddev = new((/nmons/),typeof(bw_data))
b_data_1_stddev = new((/nmons/),typeof(b_data))
bw_data_1_stddev = new((/nmons/),typeof(bw_data))
b_data_2_stddev = new((/nmons/),typeof(b_data))
bw_data_2_stddev = new((/nmons/),typeof(bw_data))

m = 0
do m = 0, nmons -1
  inds = ind((mm_b .eq. (m+1)).and.(yy_b .ge. 100.0))
  b_data_avg(m) = dim_avg_n(b_data(inds),0)
  b_data_stddev(m) = dim_stddev_n(b_data(inds),0)
  b_data_1_avg(m) = dim_avg_n(b_data_1(inds),0)
  b_data_1_stddev(m) = dim_stddev_n(b_data_1(inds),0)
  b_data_2_avg(m) = dim_avg_n(b_data_2(inds),0)
  b_data_2_stddev(m) = dim_stddev_n(b_data_2(inds),0)
  delete(inds)
  inds = ind((mm_bw .eq. (m+1)).and.(yy_bw .ge. 100.0))
  bw_data_avg(m) = dim_avg_n(bw_data(inds),0)
  bw_data_stddev(m) = dim_stddev_n(bw_data(inds),0)
  bw_data_1_avg(m) = dim_avg_n(bw_data_1(inds),0)
  bw_data_1_stddev(m) = dim_stddev_n(bw_data_1(inds),0)
  bw_data_2_avg(m) = dim_avg_n(bw_data_2(inds),0)
  bw_data_2_stddev(m) = dim_stddev_n(bw_data_2(inds),0)
  delete(inds)
end do

;;;;;;;;;;
; Combine for plotting
;;;;;;;;;;
b_data_sub = new((/3,nmons/),typeof(b_data))
bw_data_sub = new((/3,nmons/),typeof(bw_data))

b_data_sub(0,:) = b_data_avg
b_data_sub(1,:) = b_data_avg + b_data_stddev
b_data_sub(2,:) = b_data_avg - b_data_stddev

bw_data_sub(0,:) = bw_data_avg
bw_data_sub(1,:) = bw_data_avg + bw_data_stddev
bw_data_sub(2,:) = bw_data_avg - bw_data_stddev

; down
b_data_1_sub = new((/3,nmons/),typeof(b_data))
bw_data_1_sub = new((/3,nmons/),typeof(bw_data))

b_data_1_sub(0,:) = b_data_1_avg
b_data_1_sub(1,:) = b_data_1_avg + b_data_1_stddev
b_data_1_sub(2,:) = b_data_1_avg - b_data_1_stddev

bw_data_1_sub(0,:) = bw_data_1_avg
bw_data_1_sub(1,:) = bw_data_1_avg + bw_data_1_stddev
bw_data_1_sub(2,:) = bw_data_1_avg - bw_data_1_stddev

; up
b_data_2_sub = new((/3,nmons/),typeof(b_data))
bw_data_2_sub = new((/3,nmons/),typeof(bw_data))

b_data_2_sub(0,:) = b_data_2_avg
b_data_2_sub(1,:) = b_data_2_avg + b_data_2_stddev
b_data_2_sub(2,:) = b_data_2_avg - b_data_2_stddev

bw_data_2_sub(0,:) = bw_data_2_avg
bw_data_2_sub(1,:) = bw_data_2_avg + bw_data_2_stddev
bw_data_2_sub(2,:) = bw_data_2_avg - bw_data_2_stddev

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_and_terms_annual_cycle"

xdata = ispan(1, 12, 1)
label_inds = xdata
label_str = flt2string(int2flt(xdata))

if (varcode .eq. "siswnet")
   hemi_str = "Northern Hemisphere"
   ymin = -250.0
   ymax = 250.0
   yminor = 1
   x0 = 6.5
   x1 = 3.0
   x2 = x1
   x3 = x1
   y0 = 15.0
   y1 = 8.0
   y2 = 7.0
   y3 = 6.0
end if
if (varcode .eq. "silwnet")
   hemi_str = "Northern Hemisphere"
   ymin = -320.0
   ymax = 320.0
   yminor = 1
   x0 = 6.5
   x1 = 3.0
   x2 = x1
   x3 = x1
   y0 = 15.0
   y1 = 8.0
   y2 = 7.0
   y3 = 6.0
end if
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
  opts@trXMaxF                     = 12
  opts@tiXAxisString               = "Month"
  opts@tiYAxisString               = units
  opts@trYMinF                     = ymin
  opts@trYMaxF                     = ymax
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  ;opts@tmYLMinorPerMajor           = yminor
  ;opts@tmYRMinorPerMajor           = yminor
  opts@gsnYRefLine                  = 0.0

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2

; resources for b - net
  opts_b = opts
  opts_b@tiMainOn = False
  opts_b@xyMarkLineMode = "Lines"
  opts_b@xyLineColors = (/"black","black","black"/)
  opts_b@xyLineThicknessF = 3.5
  opts_b@xyDashPatterns = (/0,2,2/)
  txres_b = txres
  txres_b@txFontColor = "black"

; resources for bw - net
  opts_bw = opts
  opts_bw@tiMainOn = False
  opts_bw@xyMarkLineMode = "Lines"
  opts_bw@xyLineColors = (/"grey","grey","grey"/)
  opts_bw@xyLineThicknessF = 3.5
  opts_bw@xyDashPatterns = (/0,2,2/)
  txres_bw = txres
  txres_bw@txFontColor = "grey"

; resources for b - down
  opts_b1 = opts
  opts_b1@tiMainOn = False
  opts_b1@xyMarkLineMode = "Lines"
  opts_b1@xyLineColors = (/"firebrick","firebrick","firebrick"/)
  opts_b1@xyLineThicknessF = 3.5
  opts_b1@xyDashPatterns = (/0,2,2/)
  txres_b1 = txres
  txres_b1@txFontColor = "firebrick"

; resources for bw - down
  opts_bw1 = opts
  opts_bw1@tiMainOn = False
  opts_bw1@xyMarkLineMode = "Lines"
  opts_bw1@xyLineColors = (/"pink","pink","pink"/)
  opts_bw1@xyLineThicknessF = 3.5
  opts_bw1@xyDashPatterns = (/0,2,2/)
  txres_bw1 = txres
  txres_bw1@txFontColor = "pink"

; resources for b - up
  opts_b2 = opts
  opts_b2@tiMainOn = False
  opts_b2@xyMarkLineMode = "Lines"
  opts_b2@xyLineColors = (/"blue","blue","blue"/)
  opts_b2@xyLineThicknessF = 3.5
  opts_b2@xyDashPatterns = (/0,2,2/)
  txres_b2 = txres
  txres_b2@txFontColor = "blue"

; resources for bw - up
  opts_bw2 = opts
  opts_bw2@tiMainOn = False
  opts_bw2@xyMarkLineMode = "Lines"
  opts_bw2@xyLineColors = (/"deepskyblue","deepskyblue","deepskyblue"/)
  opts_bw2@xyLineThicknessF = 3.5
  opts_bw2@xyDashPatterns = (/0,2,2/)
  txres_bw2 = txres
  txres_bw2@txFontColor = "deepskyblue"

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
   lplot = gsn_csm_xy(wks,xdata,b_data_sub,opts_b)
   oplot1 = gsn_csm_xy(wks,xdata,bw_data_sub,opts_bw)
   overlay(lplot,oplot1)
   oplot2 = gsn_csm_xy(wks,xdata,b_data_1_sub,opts_b1)
   overlay(lplot,oplot2)
   oplot3 = gsn_csm_xy(wks,xdata,bw_data_1_sub,opts_bw1)
   overlay(lplot,oplot3)
   oplot4 = gsn_csm_xy(wks,xdata,b_data_2_sub,opts_b2)
   overlay(lplot,oplot4)
   oplot5 = gsn_csm_xy(wks,xdata,bw_data_2_sub,opts_bw2)
   overlay(lplot,oplot5)
   ;dum0 = gsn_add_text(wks,lplot,"PI-Control",x0,y0,txres_b)
   ;dum1 = gsn_add_text(wks,lplot,"CAM",x1,y1,txres_b)
   ;dum2 = gsn_add_text(wks,lplot,"WACCM",x2,y2,txres_bw)
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
