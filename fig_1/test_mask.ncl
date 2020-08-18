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
varcode = "sithick" 
hemi = "NH"
varcode2 = "siconc"
; ******************************************************
; BEGIN SCRIPT
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/"
grid_dir = "/glade/p/cesm/omwg/grids/"
nens_b = 1 ;10  ; Number of CAM ensembles

; choose years of interest
nmons = 12
; just 2001-2005
yr_st = 2001.
yr_ed = 2005.

; preallocate for files
b_files = new((/nens_b/),"string")
b_files_2 = new((/nens_b/),"string")

n = 0
do n = 0, nens_b -1
   i = n+1
   root = "r"+flt2string(int2flt(i))+"i1p1f1"
   b_files(n) = cmip_dir+"/CESM2/historical/"+root+"/SImon/"+varcode+"/gn/latest/"
   b_files_2(n) = cmip_dir+"/CESM2/historical/"+root+"/SImon/"+varcode2+"/gn/latest/"
   delete(root)
end do

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
delete([/b_file,dims/])

; get dates we want to load - note, the months are off by one based on time file was written
b_file = addfile(systemfunc("ls "+b_files(0)+varcode+"*.nc"),"r")
date_st = "1975-1" 
date_ed = "2014-12"
time_dum1 = b_file->time
time_dum2 = cd_calendar(time_dum1,0)
time_dum3 = flt2string(time_dum2(:,0))+"-"+flt2string(time_dum2(:,1))
ind1 = ind(time_dum3 .eq. date_st)
ind2 = ind(time_dum3 .eq. date_ed)
time_b = time_dum1(ind1:ind2)
time_b_all = time_dum2(ind1:ind2,:)
yy_b_all = time_b_all(:,0)
mm_b_all = time_b_all(:,1)
ntimes = dimsizes(time_b)
time_b_str = flt2string(yy_b_all)+"-"+flt2string(mm_b_all)
years_all = onedtond(yy_b_all, (/ntimes/12, 12/))
years = years_all(:,0)
nyrs = dimsizes(years)
mons = fspan(1.,12.,12)
nmons = dimsizes(mons)

delete([/b_file,time_dum1,time_dum2,time_dum3,ind1,ind2/])

;;;;;;;;;;
; Get CESM mask for where there is IceSat data
;;;;;;;;;;
; general mask for regions
fin = addfile(data_dir+"REGION_MASK_gx1v6.nc","r")
nh_mask = fin->nh_mask
nh_mask@_FillValue = -9999
nh_mask = where(ismissing(nh_mask), 0.0, nh_mask)
Lab_mask = fin->Lab_mask
Lab_mask@_FillValue = -9999
Lab_mask = where(ismissing(Lab_mask), 0.0, Lab_mask)
GIN_mask = fin->GIN_mask
GIN_mask@_FillValue = -9999
GIN_mask = where(ismissing(GIN_mask), 0.0, GIN_mask)
Bar_mask = fin->Bar_mask
Bar_mask@_FillValue = -9999
Bar_mask = where(ismissing(Bar_mask), 0.0, Bar_mask)
ArcOc_mask = fin->ArcOc_mask
ArcOc_mask@_FillValue = -9999
ArcOc_mask = where(ismissing(ArcOc_mask), 0.0, ArcOc_mask)
Sib_mask = fin->Sib_mask
Sib_mask@_FillValue = -9999
Sib_mask = where(ismissing(Sib_mask), 0.0, Sib_mask)
Beau_mask = fin->Beau_mask
Beau_mask@_FillValue = -9999
Beau_mask = where(ismissing(Beau_mask), 0.0, Beau_mask)
CArc_mask = fin->CArc_mask
CArc_mask@_FillValue = -9999
CArc_mask = where(ismissing(CArc_mask), 0.0, CArc_mask)
Bering_mask = fin->Bering_mask
Bering_mask@_FillValue = -9999
Bering_mask = where(ismissing(Bering_mask), 0.0, Bering_mask)
Okhotsk_mask = fin->Okhotsk_mask
Okhotsk_mask@_FillValue = -9999
Okhotsk_mask = where(ismissing(Okhotsk_mask), 0.0, Okhotsk_mask)
Hudson_mask = fin->Hudson_mask
Hudson_mask@_FillValue = -9999
Hudson_mask = where(ismissing(Hudson_mask), 0.0, Hudson_mask)
CAArch_mask = fin->CAArch_mask
CAArch_mask@_FillValue = -9999
CAArch_mask = where(ismissing(CAArch_mask), 0.0, CAArch_mask)

; create and combine masks - Marika region
   ;area_1 = where(((lat2d.ge.80.) .and. (lon2d.ge.243.0) .and. (lon2d.le.360.0)), 1.0, 0.0)
   ;area_2 = where(((lat2d.ge.80.) .and. (lon2d.ge.0.0) .and. (lon2d.le.100.0)), 1.0, 0.0)
   ;area_3 = where(((lat2d.ge.68.) .and. (lon2d.gt.100.0) .and. (lon2d.lt.243.0)), 1.0, 0.0)
   ;mask_IceSat= area_1+area_2+area_3

; Alex Jahn's IceSat region
   ;fin = addfile(data_dir+"CCSM4_Arctic_mask.nc","r")
   ;mask_IceSat = fin->CCSM4_Arctic_mask_no_Barents_Kara_CAA_Nares

; for iage only ; return
; Create mask for the IceAge product
   ;mask_IceSat = ArcOc_mask(0,:,:)+Bar_mask(0,:,:)

; for the annual cycle plots - Arctic Ocean only
  mask_IceSat = ArcOc_mask(0,:,:)

;delete(mask_IceSat)
;mask_IceSat = nh_mask(0,:,:)

;return

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data")
; data is in 50 or 100 year chunks

;;; b files
print("Loading b expt historical")
b_data = new((/nens_b,ntimes,nlat,nlon/),"float")
n = 0
do n = 0, nens_b -1
   print("Data for b ensemble "+(n+1))
   ; get variable of interest
   files = systemfunc("ls "+b_files(n)+varcode+"_*.nc")
   if (dimsizes(files) .eq. 1) then
      fin = addfile(files(0),"r")
      time = fin->time
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = ind(time_str .eq. date_ed)
      b_data(n,:,:,:) = fin->$varcode$(ind1:ind2,:,:)
      delete([/fin,time,time_dum,time_str/])
   else
      ; these are in 50 year chunks, so we need to load in the last two only (1950-1999; 2000-2014)
      ; first chunk
      fin = addfile(files(2),"r")
      time = fin->time
      dimt = dimsizes(time)
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = dimt - 1
      b_data(n,0:(ind2-ind1),:,:) = fin->$varcode$(ind1:ind2,:,:)
      delete([/fin,time,time_dum,time_str/])
      ; second chunk
      fin = addfile(files(3),"r")
      b_data(n,(ind2-ind1+1):ntimes-1,:,:) = fin->$varcode$(:,:,:)
      delete([/fin,ind1,ind2/])
   end if
   delete(files)
end do

; load concentration, which we need for masking
; had to do separately because the way files were divided for siconc and sivol
; differed for a single ensemble >:-|
b_data_2 = new((/nens_b,ntimes,nlat,nlon/),"float")
n = 0
do n = 0, nens_b -1
   ;print("Data for b ensemble "+(n+1))
   ; get variable of interest
   files = systemfunc("ls "+b_files_2(n)+varcode2+"_*.nc")
   if (dimsizes(files) .eq. 1) then
      fin = addfile(files(0),"r")
      time = fin->time
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = ind(time_str .eq. date_ed)
      dum_v = fin->$varcode2$(ind1:ind2,:,:)
      b_data_2(n,:,:,:) = dum_v*(conform(dum_v,mask_IceSat,(/1,2/)))
      delete([/fin,dum_v,time,time_dum,time_str/])
   else
      ; these are in 50 year chunks, so we need to load in the last two only (1950-1999; 2000-2014)
      ; first chunk
      fin = addfile(files(2),"r")
      time = fin->time
      dimt = dimsizes(time)
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = dimt - 1
      dum_v = fin->$varcode2$(ind1:ind2,:,:)
      b_data_2(n,0:(ind2-ind1),:,:) = dum_v*(conform(dum_v,mask_IceSat,(/1,2/)))
      delete([/fin,dum_v,time,time_dum,time_str/])
      ; second chunk
      fin = addfile(files(3),"r")
      dum_v = fin->$varcode2$(:,:,:)
      b_data_2(n,(ind2-ind1+1):ntimes-1,:,:) = dum_v*(conform(dum_v,mask_IceSat,(/1,2/)))
      delete([/fin,dum_v,ind1,ind2/])
   end if
   delete(files)
end do

; now mask thickness by concentration
b_data = where(b_data_2 .ge. 15.0, b_data, b_data@_FillValue)

;;;;;;;;;;
; Ensemble mean
;;;;;;;;;;
b_data_ens = dim_avg_n_Wrap(b_data,0)

;;;;;;;;;;
; Get months of interest
;;;;;;;;;;
print ("Getting monthly averages")
b_data_avg = new((/nmons,nlat,nlon/),"float")

m = 0
do m = 0, nmons -1
   mm = mons(m)
   ; average the model data
   inds = ind((mm_b_all .eq. mm).and.((yy_b_all .ge. yr_st).and.(yy_b_all .le. yr_ed)))
   b_data_avg(m,:,:) = dim_avg_n(b_data_ens(inds,:,:),0)
   delete(inds)
   delete(mm)
end do

;;;;;;;;;;
; Get seasons
;;;;;;;;;;
print ("Getting seasonal averages to match IceSat")
nseas = 2
b_data_seas = new((/nseas,nlat,nlon/),"float")

; b case
;FM average
b_data_seas(0,:,:) = dim_avg_n(b_data_avg(1:2,:,:),0)
;ON average
b_data_seas(1,:,:) = dim_avg_n(b_data_avg(9:10,:,:),0)

;;;;;;;;;;
; Assign lat/lon for plotting
;;;;;;;;;;
b_data_avg@lat2d = lat2d
b_data_avg@lon2d = lon2d

b_data_seas@lat2d = lat2d
b_data_seas@lon2d = lon2d

; for iage only  ; return
;b_data_seas(0,:,:) = where(b_data_seas(0,:,:) .lt. 1.5, 1.5, b_data_seas(0,:,:))
;b_data_seas(0,:,:) = where(b_data_seas(0,:,:) .ge. 1.75, 1.75, b_data_seas(0,:,:))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set contour levels based on varcode
main_cntrs=(/0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0/)
diff_cntrs=(/-3.0,-2.0,-1.0,-0.75,-0.5,-0.25,0.,0.25,0.5,0.75,1.0,2.0,3.0/)
seas_str = (/"Spring - FM","Autumn - ON"/)

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
  res@mpGridAndLimbOn      = False            ; Turn on lat/lon grid
  res@mpFillOn             = True             ; Enables map fill
  res@mpGridAndLimbOn      = False    ; set to true to show lat and lon lines
  res@mpLandFillColor      = "grey"
  res@mpGeophysicalLineColor      = "black"          ; Color of continents
  res@mpGeophysicalLineThicknessF = 1.5
  res@mpFillDrawOrder      = "PostDraw"     ; Fill continents after otherwise icesat has funny lines
  res@mpPerimDrawOrder     = "PostDraw"  
  res@gsnPolarLabelSpacing = 90.
  res@gsnPolarLabelFontHeightF = 0.01

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
  pnlres@gsnPanelLabelBar           = True  ; iage use False
  pnlres@lbOrientation              = "vertical"

;;;;;;;;;;;;;
; Make regular plots
;;;;;;;;;;;;;
print("Plotting seasonal data")
;set info for this plot
   ;fout = outdir+"test_mask_region"
   ; for iage only  ; return
   ;fout = outdir+"test_mask_iage"
   fout = outdir+"test_mask_ArcOc"

; open wks
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,"julies")

; set resources
  res_main                = res
  res_main@cnLevels       = main_cntrs

; do individual plots
plot = new(1, graphic)
; B data
s = 0
res_main@gsnCenterString = " "
plot(s) = gsn_csm_contour_map(wks,b_data_seas(0,:,:),res_main)

; panel plot
print("Paneling")
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,plot/])

print("Done with all plots. Hooray!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
