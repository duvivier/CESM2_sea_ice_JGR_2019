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
;varcode = "sivol"
;hemi = "SH"
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
data_dir_NH = "/glade/p/cesm/pcwg/ice/data/icesat_10prods/"
data_dir_SH = "/glade/p/cesm/pcwg/duvivier/icesat_antarctic/"
grid_dir = "/glade/p/cesm/omwg/grids/"
nens_b = 10  ; Number of CAM ensembles
nens_bw = 3  ; Number of WACCM ensembles

; choose years of interest
nmons = 12
; just 2003-2008
yr_st = 2003.
yr_ed = 2008.

; preallocate for files
b_files = new((/nens_b/),"string")
b_files_2 = new((/nens_b/),"string")
bw_files = new((/nens_bw/),"string")
bw_files_2 = new((/nens_bw/),"string")

n = 0
do n = 0, nens_b -1
   i = n+1
   root = "r"+flt2string(int2flt(i))+"i1p1f1"
   b_files(n) = cmip_dir+"/CESM2/historical/"+root+"/SImon/"+varcode+"/gn/latest/"
   b_files_2(n) = cmip_dir+"/CESM2/historical/"+root+"/SImon/"+varcode2+"/gn/latest/"
   delete(root)
end do

n = 0
do n = 0, nens_bw -1
   i = n+1
   root = "r"+flt2string(int2flt(i))+"i1p1f1"
   bw_files(n) = cmip_dir+"/CESM2-WACCM/historical/"+root+"/SImon/"+varcode+"/gn/latest/"
   bw_files_2(n) = cmip_dir+"/CESM2/historical/"+root+"/SImon/"+varcode2+"/gn/latest/"
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
kmt = b_file->KMT
; mask tarea where there is land
tarea = where(kmt .gt. 0.5, tarea, tarea@_FillValue)
; convert from cm^2 to m^2
tarea = tarea/100.0/100.0
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
if (hemi .eq. "NH")
   ; Alice tested both these and the mask Alex made is better for IceSat comparisons
   ; This file was created by Alex Jahn (2011) to match the IceSat area for CESM
   fin = addfile(data_dir+"CCSM4_Arctic_mask.nc","r")
   mask_IceSat = fin->CCSM4_Arctic_mask_no_Barents_Kara_CAA_Nares
   delete(fin)
   ; general mask for regions
   ;fin = addfile(data_dir+"REGION_MASK_gx1v6.nc","r")
   ;mask_IceSat = fin->ArcOc_mask(0,:,:)
   ;delete(fin)
end if
if (hemi .eq. "SH")
   ; general mask for regions
   fin = addfile(data_dir+"REGION_MASK_gx1v6.nc","r")
   mask_IceSat = fin->sh_mask(0,:,:)
   delete(fin)
end if

; set threshold for CESM data
if (hemi .eq. "NH")
   thresh = 15.0
end if
if (hemi .eq. "SH")
   thresh = 50.0
end if

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
b_data = where(b_data_2 .ge. thresh, b_data, b_data@_FillValue)

;;; bw files
; loop through bw files to get the data we need
bw_data = new((/nens_bw,ntimes,nlat,nlon/),"float")
n = 0
do n = 0, nens_bw -1 
   print("Data for bw ensemble "+(n+1))
   ; get variable of interest
   files = systemfunc("ls "+bw_files(n)+varcode+"_*.nc")
   if (dimsizes(files) .eq. 1) then
      fin = addfile(files(0),"r")
      time = fin->time
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = ind(time_str .eq. date_ed)
      bw_data(n,:,:,:) = fin->$varcode$(ind1:ind2,:,:)
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
      bw_data(n,0:(ind2-ind1),:,:) = fin->$varcode$(ind1:ind2,:,:)
      delete([/fin,time,time_dum,time_str/])
      ; second chunk
      fin = addfile(files(3),"r")
      bw_data(n,(ind2-ind1+1):ntimes-1,:,:) = fin->$varcode$(:,:,:)
      delete([/fin,ind1,ind2/])
   end if
   delete(files)
end do

; load concentration, which we need for masking
; had to do separately because the way files were divided for siconc and sivol
; differed for a single ensemble >:-|
bw_data_2 = new((/nens_bw,ntimes,nlat,nlon/),"float")
n = 0
do n = 0, nens_bw -1
   ;print("Data for bw ensemble "+(n+1))
   ; get variable of interest
   files = systemfunc("ls "+bw_files_2(n)+varcode2+"_*.nc")
   if (dimsizes(files) .eq. 1) then
      fin = addfile(files(0),"r")
      time = fin->time
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = ind(time_str .eq. date_ed)
      dum_v = fin->$varcode2$(ind1:ind2,:,:)
      ; mask where concentration is less than 50%
      bw_data_2(n,:,:,:) = dum_v*(conform(dum_v,mask_IceSat,(/1,2/)))
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
      bw_data_2(n,0:(ind2-ind1),:,:) = dum_v*(conform(dum_v,mask_IceSat,(/1,2/)))
      delete([/fin,dum_v,time,time_dum,time_str/])
      ; second chunk
      fin = addfile(files(3),"r")
      dum_v = fin->$varcode2$(:,:,:)
      bw_data_2(n,(ind2-ind1+1):ntimes-1,:,:) = dum_v*(conform(dum_v,mask_IceSat,(/1,2/)))
      delete([/fin,dum_v,ind1,ind2/])
   end if
   delete(files)
end do

; now mask thickness by concentration
bw_data = where(bw_data_2 .ge. thresh, bw_data, bw_data@_FillValue)

;;;;;;;;;;
; Ensemble mean
;;;;;;;;;;
b_data_ens = dim_avg_n_Wrap(b_data,0)
bw_data_ens = dim_avg_n_Wrap(bw_data,0)

;;;;;;;;;;
; Get months of interest
;;;;;;;;;;
print ("Getting monthly averages")
b_data_avg = new((/nmons,nlat,nlon/),"float")
bw_data_avg = new((/nmons,nlat,nlon/),"float")

m = 0
do m = 0, nmons -1
   mm = mons(m)
   ; average the model data
   inds = ind((mm_b_all .eq. mm).and.((yy_b_all .ge. yr_st).and.(yy_b_all .le. yr_ed)))
   b_data_avg(m,:,:) = dim_avg_n(b_data_ens(inds,:,:),0)
   bw_data_avg(m,:,:) = dim_avg_n(bw_data_ens(inds,:,:),0)
   delete(inds)
   delete(mm)
end do

;;;;;;;;;;
; Get seasons
;;;;;;;;;;
print ("Getting seasonal averages to match IceSat")
; NOTE THAT FOR 
; THE ARCTIC FM = Spring, ON = Autumn
; THE ANTARCTIC FM = Summer, ON = Spring

nseas = 2
b_data_seas = new((/nseas,nlat,nlon/),"float")
bw_data_seas = new((/nseas,nlat,nlon/),"float")

; b case
;FM average
b_data_seas(0,:,:) = dim_avg_n(b_data_avg(1:2,:,:),0)
;ON average
b_data_seas(1,:,:) = dim_avg_n(b_data_avg(9:10,:,:),0)

; bw case
;FM average
bw_data_seas(0,:,:) = dim_avg_n(bw_data_avg(1:2,:,:),0)
;ON average
bw_data_seas(1,:,:) = dim_avg_n(bw_data_avg(9:10,:,:),0)

;;;;;;;;;;
; Get IceSat data
;;;;;;;;;;
print("Loading IceSat data - 2003-2008")
; NOTE THAT FOR 
; THE ARCTIC FM = Spring, ON = Autumn
; THE ANTARCTIC FM = Summer, ON = Spring

if (hemi .eq. "NH")
  ncols = 8    
  nrows = 11118  ; 
  fin  = asciiread(data_dir_NH+"icesat_icethk_spring_regular_grid.dat",(/nrows,ncols/),"float")
  lat_fm   = fin(:,0)
  lon_fm   = fin(:,1)
  obs_fm = fin(:,7)
  obs_fm = obs_fm/100
  obs_fm@long_name = "ice thickness from icesat fir fm04"
  obs_fm@name      = "ice thickness"
  obs_fm@units     = "m"
  delete(fin)

  fin  = asciiread(data_dir_NH+"icesat_icethk_fall_regular_grid.dat",(/nrows,ncols/),"float")
  lat_on   = fin(:,0)
  lon_on   = fin(:,1)
  obs_on = fin(:,7) ; column 8 is the mean
  obs_on = obs_on/100
  obs_on@long_name = "ice thickness from icesat fir on04"
  obs_on@name      = "ice thickness"
  obs_on@units     = "m"
  delete(fin)
  tarea_obs=25000*25000 ; grid boxes are 25km by 25 km 
  tarea_obs@units="m^2"
end if

if (hemi .eq. "SH")
  ncols = 4    
  nrows = 104912  ; 
  fin  = asciiread(data_dir_SH+"summer_ICESat_gridded_mean_thickness.txt",(/nrows,ncols/),"float")
  lat_fm   = fin(:,0)
  lon_fm   = fin(:,1)
  obs_fm = fin(:,3)
  obs_fm@long_name = "ice thickness from icesat"
  obs_fm@name      = "ice thickness"
  obs_fm@units     = "m"
  delete(fin)

  fin  = asciiread(data_dir_SH+"spring_ICESat_gridded_mean_thickness.txt",(/nrows,ncols/),"float")
  lat_on   = fin(:,0)
  lon_on   = fin(:,1)
  obs_on = fin(:,3)
  obs_on@long_name = "ice thickness from icesat"
  obs_on@name      = "ice thickness"
  obs_on@units     = "m"
  delete(fin)
  tarea_obs=25000*25000 ; grid boxes are 25km by 25 km 
  tarea_obs@units="m^2"

  ; mask missing data
  obs_fm = where(obs_fm .eq. 0.0, -1.0, obs_fm)
  obs_on = where(obs_on .eq. 0.0, -1.0, obs_on)
end if

;;;;;;;;;;
; Assign lat/lon for plotting
;;;;;;;;;;
b_data_avg@lat2d = lat2d
b_data_avg@lon2d = lon2d
bw_data_avg@lat2d = lat2d
bw_data_avg@lon2d = lon2d

b_data_seas@lat2d = lat2d
b_data_seas@lon2d = lon2d
bw_data_seas@lat2d = lat2d
bw_data_seas@lon2d = lon2d

obs_fm@lat1d = lat_fm
obs_fm@lon1d = lon_fm
obs_on@lat1d = lat_on
obs_on@lon1d = lon_on

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set contour levels based on varcode
if (hemi .eq. "NH")
   main_cntrs=(/0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0/)
   seas_str = (/"Spring - FM","Autumn - ON"/)
end if
if (hemi .eq. "SH")
   main_cntrs=(/0,0.5,1.0,1.5,2.0,2.5,3.0/)
   seas_str = (/"Summer - FM","Spring - ON"/)
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
  res@gsnAddCyclic        = True          ; global data, wrap around
  
; Map resources for main dataset
  res@gsnPolar             = hemi
  if (hemi .eq. "NH")
     res@mpMinLatF            = 65.0
  end if
  if (hemi .eq. "SH")
     res@mpMaxLatF            = -55.0
  end if
  res@mpGridAndLimbOn      = False            ; Turn on lat/lon grid
  res@mpFillOn             = True             ; Enables map fill
  res@mpLandFillColor      = "grey"
  res@mpGeophysicalLineColor      = "black"          ; Color of continents
  res@mpGeophysicalLineThicknessF = 1.5
  res@mpFillDrawOrder      = "PostDraw"     ; Fill continents after otherwise icesat has funny lines
  res@mpPerimDrawOrder     = "PostDraw"  
  res@gsnPolarLabelSpacing = 45
  res@gsnPolarLabelFontHeightF = 0.02

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
   fout = outdir+"CMIP6_b_bw_icesat_"+hemi+"_masked_hist_"+varcode+"_seasonal"
   title_fin = b_data@title + " ("+b_data@units+") WACCM"

; open wks
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")

; set resources
  res_main                = res
  res_main@cnLevels       = main_cntrs

; do individual plots
plot = new(6, graphic)
; Icesat
s = 0
res_main@gsnCenterString = seas_str(s)
plot(s) = gsn_csm_contour_map(wks,obs_fm,res_main)
s = 1
res_main@gsnCenterString = seas_str(s)
plot(s) = gsn_csm_contour_map(wks,obs_on,res_main)
; BW data
s = 2
res_main@gsnCenterString = " "
plot(s) = gsn_csm_contour_map(wks,bw_data_seas(0,:,:),res_main)
s = 3
plot(s) = gsn_csm_contour_map(wks,bw_data_seas(1,:,:),res_main)
; B data
s = 4
res_main@gsnCenterString = " "
plot(s) = gsn_csm_contour_map(wks,b_data_seas(0,:,:),res_main)
s = 5
plot(s) = gsn_csm_contour_map(wks,b_data_seas(1,:,:),res_main)

; panel plot
print("Paneling")
  ; make panel plot
  ;pnlres@txString  = title_fin
  gsn_panel(wks,(/plot/),(/3,2/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,title_fin,plot/])

print("Done with all plots. Hooray!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
