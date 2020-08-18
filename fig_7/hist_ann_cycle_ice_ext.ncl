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
;hemi = "NH"     ; choose hemisphere - "NH" or "SH"

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
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/"
grid_dir = "/glade/p/cesm/omwg/grids/"
nens_b = 10  ; Number of CAM ensembles
nens_bw = 3  ; Number of WACCM ensembles

; choose years of interest
nmons = 12
; just 1979-2014
yr_st = 1979.
yr_ed = 2014.

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
   bw_files_2(n) = cmip_dir+"/CESM2-WACCM/historical/"+root+"/SImon/"+varcode2+"/gn/latest/"
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
delete(b_file)

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
nmons = 12

delete([/b_file,time_dum1,time_dum2,time_dum3,ind1,ind2/])

;;;;;;;;;;
; Get CESM mask for region
;;;;;;;;;;
print("get mask of region of interest")
;; get hemispheric mask
if (hemi .eq. "NH")then
   fin = addfile(data_dir+"REGION_MASK_gx1v6.nc","r")
   mask_in = fin->nh_mask(0,:,:)
   delete(fin)
end if
if (hemi .eq. "SH")then
   fin = addfile(data_dir+"REGION_MASK_gx1v6.nc","r")
   mask_in = fin->sh_mask(0,:,:)
   delete(fin)
end if
; reset missing value info
mask_in@_FillValue = -9999
mask_in = where(ismissing(mask_in), 0.0, mask_in)

; get hemispheric mask of areas
mask_hemi = tofloat(tarea*mask_in)

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data, Calculate hemispheric means at the same time")
; data is in 50 or 100 year chunks

;;; b files
; loop through b files to get the data we need
b_data_all = new((/nens_b,ntimes,nlat,nlon/),"float")
b_data_2_all = new((/nens_b,ntimes,nlat,nlon/),"float")
n = 0
do n = 0, nens_b -1
   print("Data for b ensemble "+(n+1))
   ; get variable of interest
   files = systemfunc("ls "+b_files(n)+varcode+"_*.nc")
   files2 = systemfunc("ls "+b_files_2(n)+varcode2+"_*.nc")
   if (dimsizes(files) .eq. 1) then
      fin = addfile(files(0),"r")
      fin2 = addfile(files2(0),"r")
      time = fin->time
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = ind(time_str .eq. date_ed)
      b_data_all(n,:,:,:) = fin->$varcode$(ind1:ind2,:,:)
      b_data_2_all(n,:,:,:) = fin2->$varcode2$(ind1:ind2,:,:)
      delete([/fin,fin2,time,time_dum,time_str/])
   else
      ; these are in 50 year chunks, so we need to load in the last two only (1950-1999; 2000-2014)
      ; first chunk
      fin = addfile(files(2),"r")
      fin2 = addfile(files2(2),"r")
      time = fin->time
      dimt = dimsizes(time)
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = dimt - 1
      b_data_all(n,0:(ind2-ind1),:,:) = fin->$varcode$(ind1:ind2,:,:)
      b_data_2_all(n,0:(ind2-ind1),:,:) = fin2->$varcode2$(ind1:ind2,:,:)
      delete([/fin,fin2,time,time_dum,time_str/])
      ; second chunk
      fin = addfile(files(3),"r")
      fin2 = addfile(files2(3),"r")
      dum_v = fin->$varcode$(:,:,:)
      b_data_all(n,(ind2-ind1+1):ntimes-1,:,:) = fin->$varcode$(:,:,:) 
      b_data_2_all(n,(ind2-ind1+1):ntimes-1,:,:) = fin2->$varcode2$(:,:,:) 
      delete([/fin,fin2,ind1,ind2/])
   end if
   delete([/files,files2/])
end do

;;; bw files
; loop through bw files to get the data we need
bw_data_all = new((/nens_bw,ntimes,nlat,nlon/),"float")
bw_data_2_all = new((/nens_bw,ntimes,nlat,nlon/),"float")
n = 0
do n = 0, nens_bw -1
   print("Data for bw ensemble "+(n+1))
   ; get variable of interest
   files = systemfunc("ls "+bw_files(n)+varcode+"_*.nc")
   files2 = systemfunc("ls "+bw_files_2(n)+varcode2+"_*.nc")
   if (dimsizes(files) .eq. 1) then
      fin = addfile(files(0),"r")
      fin2 = addfile(files2(0),"r")
      time = fin->time
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = ind(time_str .eq. date_ed)
      bw_data_all(n,:,:,:) = fin->$varcode$(ind1:ind2,:,:)
      bw_data_2_all(n,:,:,:) = fin2->$varcode2$(ind1:ind2,:,:)
      delete([/fin,fin2,time,time_dum,time_str/])
   else
      ; these are in 50 year chunks, so we need to load in the last two only (1950-1999; 2000-2014)
      ; first chunk
      fin = addfile(files(2),"r")
      fin2 = addfile(files2(2),"r")
      time = fin->time
      dimt = dimsizes(time)
      time_dum = cd_calendar(time,0)
      time_str = flt2string(time_dum(:,0))+"-"+flt2string(time_dum(:,1))
      ind1 = ind(time_str .eq. date_st)
      ind2 = dimt - 1
      bw_data_all(n,0:(ind2-ind1),:,:) = fin->$varcode$(ind1:ind2,:,:)
      bw_data_2_all(n,0:(ind2-ind1),:,:) = fin2->$varcode2$(ind1:ind2,:,:)
      delete([/fin,fin2,time,time_dum,time_str/])
      ; second chunk
      fin = addfile(files(3),"r")
      fin2 = addfile(files2(3),"r")
      dum_v = fin->$varcode$(:,:,:)
      bw_data_all(n,(ind2-ind1+1):ntimes-1,:,:) = fin->$varcode$(:,:,:) 
      bw_data_2_all(n,(ind2-ind1+1):ntimes-1,:,:) = fin2->$varcode2$(:,:,:) 
      delete([/fin,fin2,ind1,ind2/])
   end if
   delete([/files,files2/])
end do

;;;;;;;;;;
; Calculate variable over hemisphere
;;;;;;;;;;
print("calculate hemispheric extent")
; ice extent = sum of area for all grid cells with concentration > 15%
; ice extent = sum(tarea where conc>15%)

; calculate hemispheric sum 
b_data = dim_sum_n(dim_sum_n(where(b_data_2_all .ge. 15.0, 1.0, 0.0)*conform(b_data_all,mask_hemi,(/2,3/)),3),2)
bw_data = dim_sum_n(dim_sum_n(where(bw_data_2_all .ge. 15.0, 1.0, 0.0)*conform(bw_data_all,mask_hemi,(/2,3/)),3),2)

; convert to 10^6 km^2 for best y axis units
b_data = b_data/(10.^6)/(10.^6)
bw_data = bw_data/(10.^6)/(10.^6)

units = "10~S~6~N~ km~S~2~N~"

;;;;;;;;;;
; Get climatologies and ensemble means
;;;;;;;;;;
print("Getting monthly value climatologies")
b_ext_hemi_mon = new((/nens_b,nmons/),typeof(b_data))
bw_ext_hemi_mon = new((/nens_bw,nmons/),typeof(bw_data))
b_ext_hemi_avg = new((/nmons/),typeof(b_data))
bw_ext_hemi_avg = new((/nmons/),typeof(b_data))
b_ext_hemi_stddev = new((/nmons/),typeof(b_data))
bw_ext_hemi_stddev = new((/nmons/),typeof(bw_data))
b_ext_hemi_var = new((/nmons/),typeof(b_data))
bw_ext_hemi_var = new((/nmons/),typeof(bw_data))
m = 0
do m = 0, nmons -1 
   inds = ind((mm_b_all .eq. (m+1)) .and. (yy_b_all .ge. yr_st) .and. (yy_b_all .le. yr_ed))
   b_ext_hemi_mon(:,m) = dim_avg_n(b_data(:,inds),1)
   bw_ext_hemi_mon(:,m) = dim_avg_n(bw_data(:,inds),1)
   nvalues = dimsizes(inds)
   b_ext_hemi_avg(m) = avg(b_data(:,inds))
   bw_ext_hemi_avg(m) = avg(bw_data(:,inds))
   b_ext_hemi_stddev(m) = stddev(b_data(:,inds))
   bw_ext_hemi_stddev(m) = stddev(bw_data(:,inds))
   b_ext_hemi_var(m) = variance(b_data(:,inds))
   bw_ext_hemi_var(m) = variance(bw_data(:,inds))
   delete(inds)
end do

;;;;;;;;;;
; Get NSIDC extent
;;;;;;;;;;
print("Getting NSIDC data")

mm_strs = (/"01","02","03","04","05","06","07","08","09","10","11","12"/)
nsidc_all = new((/nmons,nyrs/),"float")
m = 0
do m = 0, nmons -1
   if (hemi .eq. "NH")
      fin = asciiread(data_dir+"/NSIDC_SeaIce_extent/N_"+mm_strs(m)+"_extent_v3.0.csv",(/40,4/),"float")
   end if
   if (hemi .eq. "SH")
      fin = asciiread(data_dir+"/NSIDC_SeaIce_extent/S_"+mm_strs(m)+"_extent_v3.0.csv",(/40,4/),"float")
   end if
   if (m .eq. 1)
      years_nsidc = fin(:,0)
   end if
   nsidc_all(m,:) = fin(:,2)   ; units of 10^6 km^2
end do
nsidc_all = where(nsidc_all .ge. 0, nsidc_all, nsidc_all@_FillValue)

; get monthly means over years of interest (1979-2014)
inds = ind((years_nsidc .ge. yr_st) .and. (years_nsidc .le. yr_ed))
nsidc_ext_hemi_mon = dim_avg_n(nsidc_all(:,inds),1)
nsidc_ext_hemi_stddev = dim_stddev_n(nsidc_all(:,inds),1)
nsidc_ext_hemi_var = dim_variance_n(nsidc_all(:,inds),1)
nvalues_nsidc = dimsizes(inds)

;;;;;;;;;;
; Check if means and variances are statistically different
;;;;;;;;;;
print("Checking significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
prob1 = 100.*(1. - ttest(b_ext_hemi_avg,b_ext_hemi_var,(nens_b*nvalues),bw_ext_hemi_avg,bw_ext_hemi_var,(nens_bw*nvalues),True,False))

; use f-test to see if variances are different
prob2 = 100.*(1.-ftest(b_ext_hemi_var,(nens_b*nvalues),bw_ext_hemi_var,(nens_bw*nvalues),0))

; also check if they're different than NSIDC
prob3 = 100.*(1. - ttest(b_ext_hemi_avg,b_ext_hemi_var,(nens_b*nvalues),nsidc_ext_hemi_mon,nsidc_ext_hemi_var,nvalues_nsidc,True,False))
prob4 = 100.*(1. - ttest(bw_ext_hemi_avg,bw_ext_hemi_var,(nens_bw*nvalues),nsidc_ext_hemi_mon,nsidc_ext_hemi_var,nvalues_nsidc,True,False))
prob5 = 100.*(1. - ftest(b_ext_hemi_var,(nens_b*nvalues),nsidc_ext_hemi_var,nvalues_nsidc,0))
prob6 = 100.*(1. - ftest(bw_ext_hemi_var,(nens_bw*nvalues),nsidc_ext_hemi_var,nvalues_nsidc,0))

;;;;;;;;;;
; Array of ensemble means and obs
;;;;;;;;;;
; means
all_ext_hemi_mon = new((/9,nmons/),"float")
all_ext_hemi_mon(0,:) = b_ext_hemi_avg
all_ext_hemi_mon(1,:) = bw_ext_hemi_avg
all_ext_hemi_mon(2,:) = nsidc_ext_hemi_mon
all_ext_hemi_mon(3,:) = where((prob1 .ge. 95.), b_ext_hemi_avg, all_ext_hemi_mon@_FillValue)
all_ext_hemi_mon(4,:) = where((prob1 .ge. 95.), bw_ext_hemi_avg, all_ext_hemi_mon@_FillValue)
;all_ext_hemi_mon(3,:) = where(((prob1 .ge. 95.).and.(prob3 .ge. 95.)), b_ext_hemi_avg, all_ext_hemi_mon@_FillValue)
;all_ext_hemi_mon(4,:) = where(((prob1 .ge. 95.).and.(prob4 .ge. 95.)), bw_ext_hemi_avg, all_ext_hemi_mon@_FillValue)
;all_ext_hemi_mon(5,:) = where(((prob1 .ge. 95.).and.(prob3 .lt. 95.)), b_ext_hemi_avg, all_ext_hemi_mon@_FillValue)
;all_ext_hemi_mon(6,:) = where(((prob1 .ge. 95.).and.(prob4 .lt. 95.)), bw_ext_hemi_avg, all_ext_hemi_mon@_FillValue)
;all_ext_hemi_mon(7,:) = where(((prob1 .lt. 95.).and.(prob3 .ge. 95.)), b_ext_hemi_avg, all_ext_hemi_mon@_FillValue)
;all_ext_hemi_mon(8,:) = where(((prob1 .lt. 95.).and.(prob4 .ge. 95.)), bw_ext_hemi_avg, all_ext_hemi_mon@_FillValue)

; stddev
all_ext_hemi_stddev = new((/9,nmons/),"float")
all_ext_hemi_stddev(0,:) = b_ext_hemi_stddev
all_ext_hemi_stddev(1,:) = bw_ext_hemi_stddev
all_ext_hemi_stddev(2,:) = nsidc_ext_hemi_stddev
all_ext_hemi_stddev(3,:) = where((prob2 .ge. 95.), b_ext_hemi_stddev, all_ext_hemi_stddev@_FillValue)
all_ext_hemi_stddev(4,:) = where((prob2 .ge. 95.), bw_ext_hemi_stddev, all_ext_hemi_stddev@_FillValue)
;all_ext_hemi_stddev(3,:) = where(((prob2 .ge. 95.).and.(prob5 .ge. 95.)), b_ext_hemi_stddev, all_ext_hemi_stddev@_FillValue)
;all_ext_hemi_stddev(4,:) = where(((prob2 .ge. 95.).and.(prob6 .ge. 95.)), bw_ext_hemi_stddev, all_ext_hemi_stddev@_FillValue)
;all_ext_hemi_stddev(5,:) = where(((prob2 .ge. 95.).and.(prob5 .lt. 95.)), b_ext_hemi_stddev, all_ext_hemi_stddev@_FillValue)
;all_ext_hemi_stddev(6,:) = where(((prob2 .ge. 95.).and.(prob6 .lt. 95.)), bw_ext_hemi_stddev, all_ext_hemi_stddev@_FillValue)
;all_ext_hemi_stddev(7,:) = where(((prob2 .lt. 95.).and.(prob5 .ge. 95.)), b_ext_hemi_stddev, all_ext_hemi_stddev@_FillValue)
;all_ext_hemi_stddev(8,:) = where(((prob2 .lt. 95.).and.(prob6 .ge. 95.)), bw_ext_hemi_stddev, all_ext_hemi_stddev@_FillValue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
outdir = "./"
fout0 = outdir+"CMIP6_b_bw_"+hemi+"_hist_"+vartitle+"_for_annual_cycle"

ncdf = addfile(fout0+".nc","c")
ncdf->b_ext_hemi_mon = b_ext_hemi_mon
ncdf->bw_ext_hemi_mon = bw_ext_hemi_mon
ncdf->all_ext_hemi_mon = all_ext_hemi_mon
ncdf->b_ext_hemi_stddev = b_ext_hemi_stddev
ncdf->bw_ext_hemi_stddev = bw_ext_hemi_stddev
ncdf->all_ext_hemi_stddev = all_ext_hemi_stddev

print("Output file with values")
return

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_hist_"+vartitle+"_annual_cycle"

xdata = ispan(1, 12, 1)
label_inds = xdata
label_str = flt2string(int2flt(xdata))

if (hemi .eq. "NH")
   pref = "Arctic -"
   ymin = 4.0
   ymax = 16.0
   yminor = 1
   ymin1 = 0.2
   ymax1 = 1.7
   yminor1 = 2
   x0 = 3.0 ;6.5
   x1 = 3.0
   x2 = x1
   x3 = x1
   x4 = 4.0
   y0 = 5.0 ;15.0
   y1 = 11.0
   y2 = 10.0
   y3 = 9.0
   y4 = 1.5
end if
if (hemi .eq. "SH")
   pref = "Antarctic -"
   ymin = 0.0
   ymax = 20.0
   yminor = 2
   ymin1 = 0.3
   ymax1 = 1.0
   yminor1 = 1
   x0 = 3.0 ;6.5
   x1 = 9.0
   x2 = x1
   x3 = x1
   x4 = 4.0
   y0 = 19.0
   y1 = 8.0
   y2 = 7.0
   y3 = 6.0
   y4 = 0.95
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
  opts_b@xyLineColor = "grey"
  opts_b@xyLineThicknessF = 3.5
  opts_b@xyDashPattern = 0
  txres_b = txres
  txres_b@txFontColor = "black"

; resources for bw
  opts_bw = opts
  opts_bw@tiMainOn = False
  opts_bw@xyMarkLineMode = "Lines"
  opts_bw@xyLineColor = "pink"
  opts_bw@xyLineThicknessF = 3.5
  opts_bw@xyDashPattern = 0
  txres_bw = txres
  txres_bw@txFontColor = "firebrick"

; resources for nsidc
  opts_o = opts
  opts_o@tiMainOn = False
  opts_o@xyMarkLineMode = "MarkLines"
  opts_o@xyLineColors = (/"black","firebrick","blue","black","firebrick","black","firebrick","black","firebrick"/)
  opts_o@xyLineThicknessF = 5.0
  opts_o@xyDashPatterns = (/0,0,0,0,0,0,0/)
  opts_o@xyMarkers = (/1,1,1,16,16,4,4,9,9/) 
  opts_o@xyMarkerColors = opts_o@xyLineColors
  txres_o = txres
  txres_o@txFontColor = "blue"

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
   lplot = gsn_csm_xy(wks,xdata,b_ext_hemi_mon,opts_b)
   oplot1 = gsn_csm_xy(wks,xdata,bw_ext_hemi_mon,opts_bw)
   overlay(lplot,oplot1)
   oplot2 = gsn_csm_xy(wks,xdata,all_ext_hemi_mon,opts_o)
   overlay(lplot,oplot2)
   dum0 = gsn_add_text(wks,lplot,"b) "+pref+"Historical",x0,y0,txres_b)
   dum1 = gsn_add_text(wks,lplot,"NSIDC",x1,y1,txres_o)
   dum2 = gsn_add_text(wks,lplot,"CAM6",x2,y2,txres_b)
   dum3 = gsn_add_text(wks,lplot,"WACCM6",x3,y3,txres_bw)
   plot(p) = lplot
   ;seasonal stddev
   p = 1

  opts_b@trYMinF                     = ymin1
  opts_b@trYMaxF                     = ymax1
  opts_b@tmYLMinorPerMajor           = yminor1
  opts_b@tmYRMinorPerMajor           = yminor1
  opts_bw@trYMinF                     = ymin1
  opts_bw@trYMaxF                     = ymax1
  opts_bw@tmYLMinorPerMajor           = yminor1
  opts_bw@tmYRMinorPerMajor           = yminor1

   lplot = gsn_csm_xy(wks,xdata,b_ext_hemi_stddev,opts_b)
   oplot1 = gsn_csm_xy(wks,xdata,bw_ext_hemi_stddev,opts_bw)
   overlay(lplot,oplot1)
   oplot2 = gsn_csm_xy(wks,xdata,all_ext_hemi_stddev,opts_o)
   overlay(lplot,oplot2)
   dum0 = gsn_add_text(wks,lplot,"d) "+pref+"Standard Deviation",x4,y4,txres_b)
   plot(p) = lplot

;;;;;;;;;;;;;
; Final Plots
;;;;;;;;;;;;;
  ; make panel plot
  gsn_panel(wks,(/plot/),(/2,1/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Plot 1 complete")

;;;;;;;;;;;;;;;;;;;;;; END script
end
