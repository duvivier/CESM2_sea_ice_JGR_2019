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
varcode = "aice"
; cesm directory
cesm_dir = "/glade/collections/cdg/timeseries-cmip6/"
data_dir = "/glade/p/cesm/pcwg/ice/data/NSIDC_SeaIce_extent/"
grid_dir = "/glade/p/cesm/omwg/grids/"
; b - cases
root_b = "b.e21.B1850.f09_g17.CMIP6-piControl.001"
; bw - cases
root_bw = "b.e21.BW1850.f09_g17.CMIP6-piControl.001"

; choose month
mm = 13.0   ; note: 13.0 = annual mean

; choose hemisphere
hemi = "NH"
; "NH" "SH" "both"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
b_files = systemfunc("ls "+cesm_dir+root_b+"/ice/proc/tseries/month_1/*."+varcode+".*.nc")
nfiles_b = dimsizes(b_files)
bw_files = systemfunc("ls "+cesm_dir+root_bw+"/ice/proc/tseries/month_1/*."+varcode+".*.nc")
nfiles_bw = dimsizes(bw_files)

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(b_files(0),"r")
lat2d = b_file->TLAT
lon2d = b_file->TLON
tarea = b_file->tarea
dims = dimsizes(lat2d)
nlat = dims(0)
nlon = dims(1)
delete(b_file)

; get hemispheric mask
mask_hemi = tarea

   mask_NH = where(lat2d .ge. 30.0, mask_hemi, mask_hemi@_FillValue)
   mask_SH = where(lat2d .le. -30.0, mask_hemi, mask_hemi@_FillValue)
   mask_both = mask_hemi

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
b_data_NH = new((/ntimes_b/),"float")
b_data_SH = new((/ntimes_b/),"float")
b_data_both = new((/ntimes_b/),"float")
b_time = new((/ntimes_b/),"float")

n = 0
do n = 0, nfiles_b -1
   fin = addfile(b_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
   dum_v_NH = dim_sum_n(dim_sum_n(where(dum_v .ge. 0.15, 1.0, 0.0)*conform(dum_v,mask_NH,(/1,2/)),2),1) 
   dum_v_SH = dim_sum_n(dim_sum_n(where(dum_v .ge. 0.15, 1.0, 0.0)*conform(dum_v,mask_SH,(/1,2/)),2),1)
   dum_v_both = dim_sum_n(dim_sum_n(where(dum_v .ge. 0.15, 1.0, 0.0)*conform(dum_v,mask_both,(/1,2/)),2),1)  
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
   b_data_NH(ind_st:ind_ed) = dum_v_NH
   b_data_SH(ind_st:ind_ed) = dum_v_SH
   b_data_both(ind_st:ind_ed) = dum_v_both
   b_time(ind_st:ind_ed) = dum_t
   delete([/fin,dum_t,dum_v,dum_v_NH,dum_v_SH,dum_v_both,dim_t/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data_NH = new((/ntimes_bw/),"float")
bw_data_SH = new((/ntimes_bw/),"float")
bw_data_both = new((/ntimes_bw/),"float")
bw_time = new((/ntimes_bw/),"float")

n = 0
do n = 0, nfiles_bw -1
   fin = addfile(bw_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
   dum_v_NH = dim_sum_n(dim_sum_n(where(dum_v .ge. 0.15, 1.0, 0.0)*conform(dum_v,mask_NH,(/1,2/)),2),1) 
   dum_v_SH = dim_sum_n(dim_sum_n(where(dum_v .ge. 0.15, 1.0, 0.0)*conform(dum_v,mask_SH,(/1,2/)),2),1)
   dum_v_both = dim_sum_n(dim_sum_n(where(dum_v .ge. 0.15, 1.0, 0.0)*conform(dum_v,mask_both,(/1,2/)),2),1)  
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
   bw_data_NH(ind_st:ind_ed) = dum_v_NH
   bw_data_SH(ind_st:ind_ed) = dum_v_SH
   bw_data_both(ind_st:ind_ed) = dum_v_both
   bw_time(ind_st:ind_ed) = dum_t
   delete([/fin,dum_t,dum_v,dum_v_NH,dum_v_SH,dum_v_both,dim_t/])
end do

; convert to km2 (divide by 10^6) and then convert to 10^6 (divide by 10^6 again)
b_data_NH = b_data_NH/(10.^6)/(10.^6)
bw_data_NH = bw_data_NH/(10.^6)/(10.^6)
b_data_SH = b_data_SH/(10.^6)/(10.^6)
bw_data_SH = bw_data_SH/(10.^6)/(10.^6)
b_data_both = b_data_both/(10.^6)/(10.^6)
bw_data_both = bw_data_both/(10.^6)/(10.^6)

;;;;;;;;;;
; Get annual means now
;;;;;;;;;;
; NH
b_data_NH_ann_mean = avg(b_data_NH)
bw_data_NH_ann_mean = avg(bw_data_NH)
print("NH annual mean extents over entire PI control runs")
print("CAM: "+b_data_NH_ann_mean+" *10^6 km^2")
print("WACCM: "+bw_data_NH_ann_mean+" *10^6 km^2")
b_data_NH_ann_stddev = stddev(b_data_NH)
bw_data_NH_ann_stddev = stddev(bw_data_NH)
print("NH annual stddev extents over entire PI control runs")
print("CAM: "+b_data_NH_ann_stddev+" *10^6 km^2")
print("WACCM: "+bw_data_NH_ann_stddev+" *10^6 km^2")

; SH
b_data_SH_ann_mean = avg(b_data_SH)
bw_data_SH_ann_mean = avg(bw_data_SH)
print("SH annual mean extents over entire PI control runs")
print("CAM: "+b_data_SH_ann_mean+" *10^6 km^2")
print("WACCM: "+bw_data_SH_ann_mean+" *10^6 km^2")
b_data_SH_ann_stddev = stddev(b_data_SH)
bw_data_SH_ann_stddev = stddev(bw_data_SH)
print("SH annual stddev extents over entire PI control runs")
print("CAM: "+b_data_SH_ann_stddev+" *10^6 km^2")
print("WACCM: "+bw_data_SH_ann_stddev+" *10^6 km^2")

;;;;;;;;;;;;;;;;;;;;;; END script
end
