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
;varcode = "sivol"
;hemi = "NH"     ; choose hemisphere - "NH" or "SH"

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
print("calculate hemispheric volume")
; ice vol = sum of actual ice thickness per cell times the area of that cell for all grid cells with concentration > 15%
; ice vol = sum(tarea*sivol)

; calculate hemispheric sum 
b_data = dim_sum_n(dim_sum_n(b_data_all*conform(b_data_all,mask_hemi,(/2,3/)),3),2)
bw_data = dim_sum_n(dim_sum_n(bw_data_all*conform(bw_data_all,mask_hemi,(/2,3/)),3),2)

; convert to 10^3 km^3 for best y axis units
b_data = b_data/(10.^9)/(10.^3)
bw_data = bw_data/(10.^9)/(10.^3)

units = "10~S~3~N~ km~S~3~N~"

;;;;;;;;;;
; Get PIOMAS volume
;;;;;;;;;;
if (hemi .eq. "NH")
print("Getting PIOMAS data")
; load and process PIOMAS data
data = asciiread(data_dir+"/PIOMAS.vol.daily.1979.2017.Current.v2.1.dat",(/14143,3/),"float")
; subsets
yy_pio = data(:,0)
dd_pio = data(:,1)
piomas_data = data(:,2)  ; in units of 10^3 km^3

years_pio = ispan(floattoint(min(yy_pio)),floattoint(max(yy_pio)),1)
nyrs_pio = dimsizes(years_pio)
delete(data)
end if

;;;;;;;;;;
; Get GIOMAS volume
;;;;;;;;;;
print("Getting GIOMAS data")
; load and process GIOMAS data (available for both hemispheres)
data = asciiread(data_dir+"/GIOMAS."+hemi+".vol.monthly.1979.2015.dat", (/444,3/), "float")

; subsets
yy_gio = data(:,0)
mm_gio = data(:,1)
giomas_data = data(:,2)  ; in units of 10^3 km^3

years_gio = ispan(floattoint(min(yy_gio)),floattoint(max(yy_gio)),1)
nyrs_gio = dimsizes(years_gio)
delete(data)

;;;;;;;;;;
; Get climatologies
;;;;;;;;;;
print("Getting monthly value climatologies")
b_vol_hemi_mon = new((/nens_b,nmons/),typeof(b_data))
bw_vol_hemi_mon = new((/nens_bw,nmons/),typeof(bw_data))
gio_vol_hemi_mon = new((/nmons/),typeof(giomas_data))
b_vol_hemi_avg = new((/nmons/),typeof(b_data))
bw_vol_hemi_avg = new((/nmons/),typeof(bw_data))
b_vol_hemi_stddev = new((/nmons/),typeof(b_data))
bw_vol_hemi_stddev = new((/nmons/),typeof(bw_data))
gio_vol_hemi_stddev = new((/nmons/),typeof(giomas_data))
b_vol_hemi_var = new((/nmons/),typeof(b_data))
bw_vol_hemi_var = new((/nmons/),typeof(bw_data))
gio_vol_hemi_var = new((/nmons/),typeof(giomas_data))
if (hemi .eq. "NH")
   pio_vol_hemi_mon = new((/nmons/),typeof(piomas_data))
   pio_vol_hemi_stddev = new((/nmons/),typeof(piomas_data))
   pio_vol_hemi_var = new((/nmons/),typeof(piomas_data))
   dd_st = (/1,32,60,91,121,152,182,213,244,274,305,335/)
   dd_ed = (/31,59,90,120,151,181,212,243,273,304,334,365/)
end if

m = 0
do m = 0, nmons -1 
   inds = ind((mm_b_all .eq. (m+1)) .and. (yy_b_all .ge. yr_st) .and. (yy_b_all .le. yr_ed))
   b_vol_hemi_mon(:,m) = dim_avg_n(b_data(:,inds),1)
   bw_vol_hemi_mon(:,m) = dim_avg_n(bw_data(:,inds),1)
   b_vol_hemi_avg(m) = avg(b_data(:,inds))
   bw_vol_hemi_avg(m) = avg(bw_data(:,inds))
   b_vol_hemi_stddev(m) = stddev(b_data(:,inds))
   bw_vol_hemi_stddev(m) = stddev(bw_data(:,inds))
   b_vol_hemi_var(m) = variance(b_data(:,inds))
   bw_vol_hemi_var(m) = variance(bw_data(:,inds))
   nvalues = dimsizes(inds)
   delete(inds)
   inds = ind((mm_gio .eq. (m+1)) .and. (yy_gio .ge. yr_st) .and. (yy_gio .le. yr_ed))
   gio_vol_hemi_mon(m) = dim_avg_n(giomas_data(inds),0)
   gio_vol_hemi_stddev(m) = stddev(giomas_data(inds))
   gio_vol_hemi_var(m) = variance(giomas_data(inds))
   nvalues_gio = dimsizes(inds)
   delete(inds)
   if (hemi .eq. "NH")
      inds = ind((dd_pio .ge. (dd_st(m))) .and. (dd_pio .le. (dd_ed(m))) .and. (yy_pio .ge. yr_st) .and. (yy_pio .le. yr_ed))
      pio_vol_hemi_mon(m) = dim_avg_n(piomas_data(inds),0)
      pio_vol_hemi_stddev(m) = stddev(piomas_data(inds))
      pio_vol_hemi_var(m) = variance(piomas_data(inds))
      nvalues_pio = dimsizes(inds)
      delete(inds)
   end if
end do

;;;;;;;;;;
; Check if means and variances are statistically different
;;;;;;;;;;
print("Checking significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
prob1 = 100.*(1. - ttest(b_vol_hemi_avg,b_vol_hemi_var,(nens_b*nvalues),bw_vol_hemi_avg,bw_vol_hemi_var,(nens_bw*nvalues),True,False))

; use f-test to see if variances are different
prob2 = 100.*(1.-ftest(b_vol_hemi_var,(nens_b*nvalues),bw_vol_hemi_var,(nens_bw*nvalues),0))

; also check if they're different than GIOMAS
prob3 = 100.*(1. - ttest(b_vol_hemi_avg,b_vol_hemi_var,(nens_b*nvalues),gio_vol_hemi_mon,gio_vol_hemi_var,nvalues_gio,True,False))
prob4 = 100.*(1. - ttest(bw_vol_hemi_avg,bw_vol_hemi_var,(nens_bw*nvalues),gio_vol_hemi_mon,gio_vol_hemi_var,nvalues_gio,True,False))
prob5 = 100.*(1. - ftest(b_vol_hemi_var,(nens_b*nvalues),gio_vol_hemi_var,nvalues_gio,0))
prob6 = 100.*(1. - ftest(bw_vol_hemi_var,(nens_bw*nvalues),gio_vol_hemi_var,nvalues_gio,0))

if (hemi .eq. "NH")
   ; also check if they're different than GIOMAS
   prob7 = 100.*(1. - ttest(b_vol_hemi_avg,b_vol_hemi_var,(nens_b*nvalues),pio_vol_hemi_mon,pio_vol_hemi_var,nvalues_pio,True,False))
   prob8 = 100.*(1. - ttest(bw_vol_hemi_avg,bw_vol_hemi_var,(nens_bw*nvalues),pio_vol_hemi_mon,pio_vol_hemi_var,nvalues_pio,True,False))
   prob9 = 100.*(1. - ftest(b_vol_hemi_var,(nens_b*nvalues),pio_vol_hemi_var,nvalues_pio,0))
   prob10 = 100.*(1. - ftest(bw_vol_hemi_var,(nens_bw*nvalues),pio_vol_hemi_var,nvalues_pio,0))
end if

;;;;;;;;;;
; Array of ensemble means and obs
;;;;;;;;;;
; means
all_vol_hemi_mon = new((/10,nmons/),"float")
all_vol_hemi_mon(0,:) = dim_avg_n_Wrap(b_vol_hemi_mon,0)
all_vol_hemi_mon(1,:) = dim_avg_n_Wrap(bw_vol_hemi_mon,0)
all_vol_hemi_mon(2,:) = gio_vol_hemi_mon
if (hemi .eq. "NH")
   all_vol_hemi_mon(3,:) = pio_vol_hemi_mon
end if
all_vol_hemi_mon(4,:) = where((prob1 .ge. 95.), b_vol_hemi_avg, all_vol_hemi_mon@_FillValue)
all_vol_hemi_mon(5,:) = where((prob1 .ge. 95.), bw_vol_hemi_avg, all_vol_hemi_mon@_FillValue)
;all_vol_hemi_mon(4,:) = where(((prob1 .ge. 95.).and.(prob3 .ge. 95.)), b_vol_hemi_avg, all_vol_hemi_mon@_FillValue)
;all_vol_hemi_mon(5,:) = where(((prob1 .ge. 95.).and.(prob4 .ge. 95.)), bw_vol_hemi_avg, all_vol_hemi_mon@_FillValue)
;all_vol_hemi_mon(6,:) = where(((prob1 .ge. 95.).and.(prob3 .lt. 95.)), b_vol_hemi_avg, all_vol_hemi_mon@_FillValue)
;all_vol_hemi_mon(7,:) = where(((prob1 .ge. 95.).and.(prob4 .lt. 95.)), bw_vol_hemi_avg, all_vol_hemi_mon@_FillValue)
;all_vol_hemi_mon(8,:) = where(((prob1 .lt. 95.).and.(prob3 .ge. 95.)), b_vol_hemi_avg, all_vol_hemi_mon@_FillValue)
;all_vol_hemi_mon(9,:) = where(((prob1 .lt. 95.).and.(prob4 .ge. 95.)), bw_vol_hemi_avg, all_vol_hemi_mon@_FillValue)

; stddev
all_vol_hemi_stddev = new((/10,nmons/),"float")
all_vol_hemi_stddev(0,:) = b_vol_hemi_stddev
all_vol_hemi_stddev(1,:) = bw_vol_hemi_stddev
all_vol_hemi_stddev(2,:) = gio_vol_hemi_stddev
if (hemi .eq. "NH")
   all_vol_hemi_stddev(3,:) = pio_vol_hemi_stddev
end if
all_vol_hemi_stddev(4,:) = where((prob1 .ge. 95.), b_vol_hemi_stddev, all_vol_hemi_stddev@_FillValue)
all_vol_hemi_stddev(5,:) = where((prob1 .ge. 95.), bw_vol_hemi_stddev, all_vol_hemi_stddev@_FillValue)
;all_vol_hemi_stddev(4,:) = where(((prob1 .ge. 95.).and.(prob3 .ge. 95.)), b_vol_hemi_stddev, all_vol_hemi_stddev@_FillValue)
;all_vol_hemi_stddev(5,:) = where(((prob1 .ge. 95.).and.(prob4 .ge. 95.)), bw_vol_hemi_stddev, all_vol_hemi_stddev@_FillValue)
;all_vol_hemi_stddev(6,:) = where(((prob1 .ge. 95.).and.(prob3 .lt. 95.)), b_vol_hemi_stddev, all_vol_hemi_stddev@_FillValue)
;all_vol_hemi_stddev(7,:) = where(((prob1 .ge. 95.).and.(prob4 .lt. 95.)), bw_vol_hemi_stddev, all_vol_hemi_stddev@_FillValue)
;all_vol_hemi_stddev(8,:) = where(((prob1 .lt. 95.).and.(prob3 .ge. 95.)), b_vol_hemi_stddev, all_vol_hemi_stddev@_FillValue)
;all_vol_hemi_stddev(9,:) = where(((prob1 .lt. 95.).and.(prob4 .ge. 95.)), bw_vol_hemi_stddev, all_vol_hemi_stddev@_FillValue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
outdir = "./"
fout0 = outdir+"CMIP6_b_bw_"+hemi+"_hist_"+vartitle+"_for_annual_cycle"

ncdf = addfile(fout0+".nc","c")
ncdf->b_vol_hemi_mon = b_vol_hemi_mon
ncdf->bw_vol_hemi_mon = bw_vol_hemi_mon
ncdf->all_vol_hemi_mon = all_vol_hemi_mon
ncdf->b_vol_hemi_stddev = b_vol_hemi_stddev
ncdf->bw_vol_hemi_stddev = bw_vol_hemi_stddev
ncdf->all_vol_hemi_stddev = all_vol_hemi_stddev

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
   ymin = 5.0
   ymax = 30.0
   yminor = 4
   ymin1 = 2.25
   ymax1 = 6.0
   yminor1 = 3
   x0 = 9.0
   x1 = 3.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 4.0
   y0 = 29.0
   y1 = 12.5
   y2 = 11.
   y3 = 9.5
   y4 = 14.0
   y5 = 5.75
end if
if (hemi .eq. "SH")
   pref = "Antarctic -"
   ymin = 2.0
   ymax = 21.0
   yminor = 3
   ymin1 = 0.4
   ymax1 = 1.4
   yminor1 = 3
   x0 = 3.0
   x1 = 9.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 4.0
   y0 = 20.0
   y1 = 8.0
   y2 = 7.0
   y3 = 6.0
   y4 = 9.0
   y5 = 1.25 
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
  opts_o@xyLineColors = (/"black","firebrick","blue","deepskyblue","black","firebrick","black","firebrick","black","firebrick"/)
  opts_o@xyLineThicknessF = 5.0
  opts_o@xyDashPattern = 0
  opts_o@xyMarkers = (/1,1,1,1,16,16,4,4,9,9/)
  opts_o@xyMarkerColors = opts_o@xyLineColors
  txres_o = txres
  txres_o@txFontColor = "blue"
  txres_o0 = txres
  txres_o0@txFontColor = "deepskyblue"

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

   ; do actual plotting
   p = 0
   lplot = gsn_csm_xy(wks,xdata,b_vol_hemi_mon,opts_b)
   oplot1 = gsn_csm_xy(wks,xdata,bw_vol_hemi_mon,opts_bw)
   overlay(lplot,oplot1)
   oplot2 = gsn_csm_xy(wks,xdata,all_vol_hemi_mon,opts_o)
   overlay(lplot,oplot2)
   dum0 = gsn_add_text(wks,lplot,"b) "+pref+"Historical",x0,y0,txres_b)
   dum1 = gsn_add_text(wks,lplot,"GIOMAS",x1,y1,txres_o)
   dum2 = gsn_add_text(wks,lplot,"CAM6",x2,y2,txres_b)
   dum3 = gsn_add_text(wks,lplot,"WACCM6",x3,y3,txres_bw)
   if (hemi .eq. "NH")
      dum4 = gsn_add_text(wks,lplot,"PIOMAS",x4,y4,txres_o0)
   end if
   plot(p) = lplot
   ;seasonal stddev
   p = 1

  opts_o@trYMinF                     = ymin1
  opts_o@trYMaxF                     = ymax1
  opts_o@tmYLMinorPerMajor           = yminor1
  opts_o@tmYRMinorPerMajor           = yminor1

   lplot = gsn_csm_xy(wks,xdata,all_vol_hemi_stddev,opts_o)
   dum5 = gsn_add_text(wks,lplot,"d) "+pref+"Standard Deviation",x5,y5,txres_b)
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
