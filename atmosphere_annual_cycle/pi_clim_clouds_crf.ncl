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
;varcode = "swcrf"
;hemi = "NH"     ; choose hemisphere - "NH" or "SH"
; ******************************************************
; BEGIN SCRIPT

if (varcode .eq. "swcrf")
   varcode1 = "rsds"
   varcode2 = "rsdscs"
end if
if (varcode .eq. "lwcrf")
   varcode1 = "rlds"
   varcode2 = "rldscs"
end if
varcode_i = "ICEFRAC"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/NSIDC_SeaIce_extent/"
cmip_dir_i = "/glade/collections/cdg/timeseries-cmip6/"
grid_dir = "/glade/p/cesm/omwg/grids/"

b_files1 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/Amon/"+varcode1+"/gn/latest/*.nc")
b_files2 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/Amon/"+varcode2+"/gn/latest/*.nc")
b_files_i = systemfunc("ls "+cmip_dir_i+"/b.e21.B1850.f09_g17.CMIP6-piControl.001/atm/proc/tseries/month_1/*"+varcode_i+"*.nc")
nfiles_b = dimsizes(b_files1)
bw_files1 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/Amon/"+varcode1+"/gn/latest/*.nc")
bw_files2 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/Amon/"+varcode2+"/gn/latest/*.nc")
bw_files_i = systemfunc("ls "+cmip_dir_i+"/b.e21.BW1850.f09_g17.CMIP6-piControl.001/atm/proc/tseries/month_1/*"+varcode_i+"*.nc")
nfiles_bw = dimsizes(bw_files1)

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(b_files1(0),"r")
dummy2d = b_file->$varcode1$(0,:,:)
lat1d = b_file->lat
lon1d = b_file->lon
lat2d = conform(dummy2d,lat1d,0)
copy_VarMeta(dummy2d,lat2d)
lat2d@units = "degrees"
lon2d = conform(dummy2d,lon1d,1)
delete([/b_file,dummy2d/])

; get hemispheric mask
mask_hemi = tofloat(lat2d)

if (hemi .eq. "NH")then
   mask_hemi = where(lat2d .ge. 70.0, 1.0, mask_hemi@_FillValue)
end if
if (hemi .eq. "SH")then
   mask_hemi = where(lat2d .le. -50.0, 1.0, mask_hemi@_FillValue)
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
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   ; get variable of interest
   fin1 = addfile(b_files1(n),"r")
   dum_1 = fin1->$varcode1$
   fin2 = addfile(b_files2(n),"r")
   dum_2 = fin2->$varcode2$
   dum_v = dum_1 - dum_2
   ; get ice fraction on atm grid (0-1)
   fin_i = addfile(b_files_i(n),"r")
   dum_i = fin_i->$varcode_i$
   ; get hemispheric avg where there is sea ice
   dum_v_hemi = dim_avg_n(dim_avg_n(where(dum_i .ge. 0.15, dum_v, dum_v@_FillValue)*conform(dum_v,mask_hemi,(/1,2/)),2),1)
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
   b_data(ind_st:ind_ed) = dum_v_hemi
   b_time(ind_st:ind_ed) = dum_t
   delete([/fin1,fin2,fin_i,dum_t,dum_i,dum_1,dum_2,dum_v,dum_v_hemi,dim_t/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data = new((/ntimes_bw/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   ; get variable of interest
   fin1 = addfile(bw_files1(n),"r")
   dum_1 = fin1->$varcode1$
   fin2 = addfile(bw_files2(n),"r")
   dum_2 = fin2->$varcode2$
   dum_v = dum_1 - dum_2
   ; get ice fraction on atm grid (0-1)
   fin_i = addfile(bw_files_i(n),"r")
   dum_i = fin_i->$varcode_i$
   ; get hemispheric avg where there is sea ice
   dum_v_hemi = dim_avg_n(dim_avg_n(where(dum_i .ge. 0.15, dum_v, dum_v@_FillValue)*conform(dum_v,mask_hemi,(/1,2/)),2),1)
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
   bw_data(ind_st:ind_ed) = dum_v_hemi
   bw_time(ind_st:ind_ed) = dum_t
   delete([/fin1,fin2,fin_i,dum_t,dum_i,dum_1,dum_2,dum_v,dum_v_hemi,dim_t/])
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
b_data_stddev = new((/nmons/),typeof(b_data))
bw_data_stddev = new((/nmons/),typeof(bw_data))

m = 0
do m = 0, nmons -1
  inds = ind((mm_b .eq. (m+1)).and.(yy_b .ge. 100.0))
  b_data_avg(m) = dim_avg_n(b_data(inds),0)
  b_data_stddev(m) = dim_stddev_n(b_data(inds),0)
  delete(inds)
  inds = ind((mm_bw .eq. (m+1)).and.(yy_bw .ge. 100.0))
  bw_data_avg(m) = dim_avg_n(bw_data(inds),0)
  bw_data_stddev(m) = dim_stddev_n(bw_data(inds),0)
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_annual_cycle_seaice_only"

xdata = ispan(1, 12, 1)
label_inds = xdata
label_str = flt2string(int2flt(xdata))

if (hemi .eq. "NH")
   hemi_str = "Northern Hemisphere"
   ymin = 4.0
   ymax = 16.0
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
if (hemi .eq. "SH")
   hemi_str = "Southern Hemisphere"
   ymin = 0.0
   ymax = 20.0
   yminor = 2
   x0 = 6.5
   x1 = 9.0
   x2 = x1
   x3 = x1
   y0 = 19.0
   y1 = 8.0
   y2 = 7.0
   y3 = 6.0
end if

if (varcode .eq. "clt")
  units = "%"
else
  units = "W m~S~-2~N~"
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
  opts@trXMaxF                     = 12
  opts@tiXAxisString               = "Month"
  opts@tiYAxisString               = units
  ;opts@trYMinF                     = ymin
  ;opts@trYMaxF                     = ymax
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  ;opts@tmYLMinorPerMajor           = yminor
  ;opts@tmYRMinorPerMajor           = yminor

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
  opts_b@xyMarkLineMode = "Lines"
;  opts_b@xyLineColor = "black"
  opts_b@xyLineColors = (/"black","grey","grey"/)
  opts_b@xyLineThicknessF = 3.5
  opts_b@xyDashPatterns = (/0,2,2/)
  txres_b = txres
  txres_b@txFontColor = "black"

; resources for bw
  opts_bw = opts
  opts_bw@tiMainOn = False
  opts_bw@xyMarkLineMode = "Lines"
;  opts_bw@xyLineColor = "firebrick"
  opts_bw@xyLineColors = (/"firebrick","pink","pink"/)
  opts_bw@xyLineThicknessF = 3.5
  opts_bw@xyDashPatterns = (/0,2,2/)
  txres_bw = txres
  txres_bw@txFontColor = "firebrick"

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
   dum0 = gsn_add_text(wks,lplot,"PI-Control",x0,y0,txres_b)
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
