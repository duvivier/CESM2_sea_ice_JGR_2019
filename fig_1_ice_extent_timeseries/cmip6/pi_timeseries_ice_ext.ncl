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
;varcode = "siconc"
;hemi = "SH"     ; choose hemisphere - "NH" or "SH"

mm = 13.0       ; Choose month 1-12 or 13, which is annual mean
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/NSIDC_SeaIce_extent/"
grid_dir = "/glade/p/cesm/omwg/grids/"

b_files = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode+"/gn/latest/*.nc")
nfiles_b = dimsizes(b_files)
bw_files = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode+"/gn/latest/*.nc")
nfiles_bw = dimsizes(bw_files)

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
mask_hemi = tofloat(tarea)

if (hemi .eq. "NH")then
   mask_hemi = where(lat2d .ge. 30.0, mask_hemi, mask_hemi@_FillValue)
end if
if (hemi .eq. "SH")then
   mask_hemi = where(lat2d .le. -30.0, mask_hemi, mask_hemi@_FillValue)
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
; ice ext = grid cells with concentration > 15%, area of these cells (m^2)
; sum over just hemi cells

;;; b files
print("Loading b expt PI control")
b_data = new((/ntimes_b/),"float")
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   fin = addfile(b_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
   dum_v_hemi = dim_sum_n(dim_sum_n(where(dum_v .ge. 15.0, 1.0, 0.0)*conform(dum_v,mask_hemi,(/1,2/)),2),1) 
   ; deal with time index
   dum_t = fin->time
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
   delete([/fin,dum_t,dum_v,dum_v_hemi,dim_t/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data = new((/ntimes_bw/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   fin = addfile(bw_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
   dum_v_hemi = dim_sum_n(dim_sum_n(where(dum_v .ge. 15.0, 1.0, 0.0)*conform(dum_v,mask_hemi,(/1,2/)),2),1) 
   ; deal with time index
   dum_t = fin->time
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
   delete([/fin,dum_t,dum_v,dum_v_hemi,dim_t/])
end do

; convert to km2
b_data = b_data/(10.^6)
bw_data = bw_data/(10.^6)

; convert to 10^6 km^2 for better y axis units 
b_data = b_data/(10.^6)
bw_data = bw_data/(10.^6)

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
; Get just time of interest
;;;;;;;;;;
if (mm .le. 9)then
   mm_str = "0"+flt2string(mm)
else
   mm_str = flt2string(mm)
end if
; if we want the annual mean, set this
if (mm .eq. 13.0)then
   mm_str = "annual"
end if
print ("Subsetting for just "+mm_str)

; for just one month
if (mm .lt. 13.0)
   inds = ind(mm_b .eq. mm)
   b_data_sub = b_data(inds)
   yy_b_sub = yy_b(inds)
   mm_b_sub = mm_b(inds)
   delete(inds)

   inds = ind(mm_bw .eq. mm)
   bw_data_sub = bw_data(inds)
   yy_bw_sub = yy_bw(inds)
   mm_bw_sub = mm_bw(inds)
   delete(inds)
end if

; for annual mean
if (mm .eq. 13.0)
   inds = ind(mm_b .eq. 1.0)
   yy_b_sub = yy_b(inds)
   b_data_tmp = onedtond(b_data,(/nyrs_b,nmons/))
   b_data_sub = dim_avg_n_Wrap(b_data_tmp,1)
   delete([/inds,b_data_tmp/])

   inds = ind(mm_bw .eq. 1.0)
   yy_bw_sub = yy_bw(inds)
   bw_data_tmp = onedtond(bw_data,(/nyrs_bw,nmons/))
   bw_data_sub = dim_avg_n_Wrap(bw_data_tmp,1)
   delete([/inds,bw_data_tmp/])
end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_ice_ext"

label_inds = (yy_b_sub(1:nyrs_b-1:100))-2
label_minor_inds = label_inds
label_str = flt2string(label_inds)

if ((hemi .eq. "NH") .and. (mm .eq. 13.0))
   hemi_str = "Northern Hemisphere"
   ymin = 9.0
   ymax = 14.0
   yminor = 1
   x0 = 550
   x1 = 800
   x2 = x1
   y0 = 13.75
   y1 = 12.5
   y2 = 12.25
end if
if ((hemi .eq. "SH") .and. (mm .eq. 13.0))
   hemi_str = "Southern Hemisphere"
   ymin = 9.0
   ymax = 16.0
   yminor = 1
   x0 = 550
   x1 = 800
   x2 = x1
   y0 = 15.75
   y1 = 15.25
   y2 = 15.0
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
  opts@tmXBMinorValues             = label_minor_inds
  opts@trXMinF                     = 100
  opts@trXMaxF                     = 1200
  opts@tiXAxisString               = "Year"
  opts@tiYAxisString               = "10~S~6~N~ km~S~2~N~"
  opts@trYMinF                     = ymin
  opts@trYMaxF                     = ymax
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
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

; resources for b
  opts_b = opts
  opts_b@tiMainOn = False
  opts_b@xyMarkLineMode = "Lines"
  opts_b@xyLineColor = "black"
  opts_b@xyLineThicknessF = 3.5
  opts_b@xyDashPattern = 0
  txres_b = txres
  txres_b@txFontColor = "black"

; resources for bw
  opts_bw = opts
  opts_bw@tiMainOn = False
  opts_bw@xyMarkLineMode = "Lines"
  opts_bw@xyLineColor = "firebrick"
  opts_bw@xyLineThicknessF = 3.5
  opts_bw@xyDashPattern = 0
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
   lplot = gsn_csm_xy(wks,yy_b_sub,b_data_sub,opts_b)
   oplot1 = gsn_csm_xy(wks,yy_bw_sub,bw_data_sub,opts_bw)
   overlay(lplot,oplot1)
   dum1 = gsn_add_text(wks,lplot,"PI-Control",x0,y0,txres_b)
   ;dum2 = gsn_add_text(wks,lplot,"CAM",x1,y1,txres_b)
   ;dum3 = gsn_add_text(wks,lplot,"WACCM",x2,y2,txres_bw)
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
