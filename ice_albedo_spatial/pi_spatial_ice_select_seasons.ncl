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
;varcode = "siconc"
;hemi = "NH"
; ******************************************************
; BEGIN SCRIPT
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
   mask_hemi = where(lat2d .ge. 60.0, mask_hemi, mask_hemi@_FillValue)
end if
if (hemi .eq. "SH")then
   mask_hemi = where(lat2d .le. -50.0, mask_hemi, mask_hemi@_FillValue)
end if

mons = fspan(1.,12.,12)
nmons = dimsizes(mons)
nyrs_b = 1200
nyrs_bw = 499
ntimes_b = nmons*nyrs_b
ntimes_bw = nmons*nyrs_bw

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data")
; data is in 50 or 100 year chunks

;;; b files
print("Loading b expt PI control")
b_data = new((/ntimes_b,nlat,nlon/),"float")
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   fin = addfile(b_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
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
   b_data(ind_st:ind_ed,:,:) = dum_v
   b_time(ind_st:ind_ed) = dum_t
   delete([/fin,dum_t,dum_v,dim_t/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data = new((/ntimes_bw,nlat,nlon/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   fin = addfile(bw_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
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
   bw_data(ind_st:ind_ed,:,:) = dum_v
   bw_time(ind_st:ind_ed) = dum_t
   delete([/fin,dum_t,dum_v,dim_t/])
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
; Get monthly of interest
;;;;;;;;;;
print ("Getting monthly averages")
b_data_avg = new((/nmons,nlat,nlon/),"float")
bw_data_avg = new((/nmons,nlat,nlon/),"float")

m = 0
do m = 0, nmons -1
   mm = mons(m)
   ; get data
   inds = ind((mm_b .eq. mm).and.(yy_b .ge. 100.0))
   b_data_avg(m,:,:) = dim_avg_n(b_data(inds,:,:),0)
   delete(inds)
   inds = ind((mm_bw .eq. mm).and.(yy_bw .ge. 100.0))
   bw_data_avg(m,:,:) = dim_avg_n(bw_data(inds,:,:),0)
   delete(inds)   
   delete(mm)
end do

;;;;;;;;;;
; Get seasons
;;;;;;;;;;
print ("Getting seasonal averages")
nseas = 4
b_data_seas = new((/nseas,nlat,nlon/),"float")
bw_data_seas = new((/nseas,nlat,nlon/),"float")

; b case
;JFM average
b_data_seas(0,:,:) = dim_avg_n(b_data_avg(0:2,:,:),0)
;AMJ average
b_data_seas(1,:,:) = dim_avg_n(b_data_avg(3:5,:,:),0)
;JAS average
b_data_seas(2,:,:) = dim_avg_n(b_data_avg(6:8,:,:),0)
;OND average
b_data_seas(3,:,:) = dim_avg_n(b_data_avg(9:11,:,:),0)

; bw case
;JFM average
bw_data_seas(0,:,:) = dim_avg_n(bw_data_avg(0:2,:,:),0)
;AMJ average
bw_data_seas(1,:,:) = dim_avg_n(bw_data_avg(3:5,:,:),0)
;JAS average
bw_data_seas(2,:,:) = dim_avg_n(bw_data_avg(6:8,:,:),0)
;OND average
bw_data_seas(3,:,:) = dim_avg_n(bw_data_avg(9:11,:,:),0)

;;;;;;;;;;
; Get differences
;;;;;;;;;;
if (varcode .eq. "siconc")
   b_data_avg = where(b_data_avg .ge. 5.0, b_data_avg, b_data_avg@_FillValue)
   bw_data_avg = where(bw_data_avg .ge. 5.0, bw_data_avg, bw_data_avg@_FillValue)
   b_data_seas = where(b_data_seas .ge. 5.0, b_data_seas, b_data_seas@_FillValue)
   bw_data_seas = where(bw_data_seas .ge. 5.0, bw_data_seas, bw_data_seas@_FillValue)
end if
if (varcode .eq. "sivol")
   b_data_avg = where(b_data_avg .ge. 0.1, b_data_avg, b_data_avg@_FillValue)
   bw_data_avg = where(bw_data_avg .ge. 0.1, bw_data_avg, bw_data_avg@_FillValue)
   b_data_seas = where(b_data_seas .ge. 0.1, b_data_seas, b_data_seas@_FillValue)
   bw_data_seas = where(bw_data_seas .ge. 0.1, bw_data_seas, bw_data_seas@_FillValue)
end if

; calculate difference
diff_avg = b_data_avg - bw_data_avg

; calculate difference
diff_seas = b_data_seas - bw_data_seas

;;;;;;;;;;
; Assign lat/lon for plotting
;;;;;;;;;;
b_data_avg@lat2d = lat2d
b_data_avg@lon2d = lon2d
bw_data_avg@lat2d = lat2d
bw_data_avg@lon2d = lon2d
diff_avg@lat2d = lat2d
diff_avg@lon2d = lon2d

b_data_seas@lat2d = lat2d
b_data_seas@lon2d = lon2d
bw_data_seas@lat2d = lat2d
bw_data_seas@lon2d = lon2d
diff_seas@lat2d = lat2d
diff_seas@lon2d = lon2d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set contour levels based on varcode
if (varcode .eq. "siconc")
   main_cntrs=(/5,10,15,20,30,40,50, 60, 70, 80,85,90,95, 99/)
   diff_cntrs=(/-25.,-20.,-15.,-10.,-5.,0,5.,10.,15.,20.,25./)
end if
if (varcode .eq. "sivol")
   main_cntrs=(/0,0.3,0.6,1.0,1.4,1.9,2.4,3.0,3.6/)
   diff_cntrs=(/-3.0,-2.0,-1.0,-0.75,-0.5,-0.25,0.,0.25,0.5,0.75,1.0,2.0,3.0/)
end if
mons_str = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct","Nov","Dec"/)
if (hemi .eq. "NH")
   seas_str = (/"Winter - JFM","Summer - JAS"/)
end if
if (hemi .eq. "SH")
   seas_str = (/"Summer - JFM","Winter - JAS"/)
end if

; set outdir and fout names
   outdir = "./"

; set general resources
  res                     = True          ; plot mods desired
  res@gsnDraw             = False         ; don't draw the plot yet
  res@gsnFrame            = False         ; or the frame
  res@cnLinesOn           = False         ; turn off contour lines
  res@txFontHeightF       = 0.030          ; set subtitle text font
  res@gsnLeftString       = ""
  res@gsnRightString      = ""
  res@gsnCenterString     = ""
  
; Map resources for main dataset
  res@gsnPolar             = hemi
  if (hemi .eq. "NH")
     res@mpMinLatF            = 50.0
  end if
  if (hemi .eq. "SH")
     res@mpMaxLatF            = -50.0
  end if
  res@mpGridAndLimbOn      = False            ; Turn on lat/lon grid
  res@mpFillOn             = True             ; Enables map fill
  res@mpGridAndLimbOn      = False    ; set to true to show lat and lon lines
  res@mpLandFillColor      = "grey"
  res@mpGeophysicalLineColor      = "black"          ; Color of continents
  res@mpGeophysicalLineThicknessF = 1.5

  ; contour resources for everything
  res@cnFillOn            = True          ; turn on color
  res@cnFillMode          = "AreaFill"
  res@cnLinesOn           = False
  res@lbLabelBarOn        = False         ; turn off label bars
  res@cnInfoLabelOn       = False         ; turn off contour info labels
  res@cnLineLabelsOn      = False 
  res@cnMissingValFillColor = "white"
  res@cnLevelSelectionMode    = "ExplicitLevels"

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True	       ; fill the page
  pnlres@gsnPanelLabelBar           = True
  pnlres@lbOrientation              = "vertical"

;;;;;;;;;;;;;
; Make regular plots
;;;;;;;;;;;;;
print("Plotting seasonal data")
;set info for this plot
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_seasonal"
   title_fin = b_data@title + " ("+b_data@units+") WACCM"

; open wks
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,"julies")

; set resources
  res_main                = res
  res_main@cnLevels       = main_cntrs

; do individual plots
plot = new(2, graphic)
; Winter
s = 0
res_main@gsnCenterString = seas_str(s)
plot(s) = gsn_csm_contour_map(wks,bw_data_seas(0,:,:),res_main)
; Summer
s = 1
res_main@gsnCenterString = seas_str(s)
plot(s) = gsn_csm_contour_map(wks,bw_data_seas(2,:,:),res_main)

; panel plot
print("Paneling")
  ; make panel plot
  ;pnlres@txString  = title_fin
  gsn_panel(wks,(/plot/),(/1,2/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,title_fin,plot/])

;;;;;;;;;;;;;
; Make difference plots
;;;;;;;;;;;;;
print("Plotting seasonal difference data")
;set info for this plot
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_diff_seasonal"
   title_fin = b_data@title + " ("+b_data@units+") difference: CAM-WACCM"

; open wks
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,"GMT_polar")

; set resources
  res_diff                = res
  res_diff@cnLevels       = diff_cntrs

; do individual plots
plot = new(2, graphic)
; Winter
s = 0
res_diff@gsnCenterString = seas_str(s)
plot(s) = gsn_csm_contour_map(wks,diff_seas(0,:,:),res_diff)
; Summer
s = 1
res_diff@gsnCenterString = seas_str(s)
plot(s) = gsn_csm_contour_map(wks,diff_seas(2,:,:),res_diff)

; panel plot
print("Paneling")
  ; make panel plot
  ;pnlres@txString  = title_fin
  gsn_panel(wks,(/plot/),(/1,4/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,title_fin,plot/])

print("Done with all plots. Hooray!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
