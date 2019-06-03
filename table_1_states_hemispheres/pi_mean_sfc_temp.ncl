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
varcode = "ts"
; ******************************************************
varcode = "ts"
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/NSIDC_SeaIce_extent/"
grid_dir = "/glade/p/cesm/omwg/grids/"

b_files = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/Amon/"+varcode+"/gn/latest/*.nc")
nfiles_b = dimsizes(b_files)
bw_files = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/Amon/"+varcode+"/gn/latest/*.nc")
nfiles_bw = dimsizes(bw_files)

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(b_files(0),"r")
dummy2d = b_file->$varcode$(0,:,:)
lat1d = b_file->lat
lon1d = b_file->lon
lat2d = conform(dummy2d,lat1d,0)
copy_VarMeta(dummy2d,lat2d)
lat2d@units = "degrees"
lon2d = conform(dummy2d,lon1d,1)
delete([/b_file,dummy2d/])

nmons = 12
nyrs_b = 1200
nyrs_bw = 499
ntimes_b = nmons*nyrs_b
ntimes_bw = nmons*nyrs_bw

; get hemispheric mask
mask_hemi = tofloat(lat2d)
mask_NH = where(lat2d .ge. 70.0, 1.0, mask_hemi@_FillValue)
mask_NH@_FillValue = mask_hemi@_FillValue
mask_SH = where(lat2d .le. -50.0, 1.0, mask_hemi@_FillValue)
mask_SH@_FillValue = mask_hemi@_FillValue

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data, Calculate hemispheric means at the same time")
; data is in 50 or 100 year chunks
; avg over just hemi cells

;;; b files
print("Loading b expt PI control")
b_data = new((/ntimes_b/),"float")
b_data_NH = new((/ntimes_b/),"float")
b_data_SH = new((/ntimes_b/),"float")
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   fin = addfile(b_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
   dum_v_NH = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_NH,(/1,2/)),2),1) 
   dum_v_SH = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_SH,(/1,2/)),2),1)
   dum_v_both = dim_avg_n(dim_avg_n(dum_v,2),1)
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
   b_data(ind_st:ind_ed) = dum_v_both
   b_data_NH(ind_st:ind_ed) = dum_v_NH
   b_data_SH(ind_st:ind_ed) = dum_v_SH
   b_time(ind_st:ind_ed) = dum_t
   delete([/fin,dum_t,dum_v,dum_v_both,dum_v_NH,dum_v_SH,dim_t/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data = new((/ntimes_bw/),"float")
bw_data_NH = new((/ntimes_bw/),"float")
bw_data_SH = new((/ntimes_bw/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   fin = addfile(bw_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
   dum_v_NH = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_NH,(/1,2/)),2),1) 
   dum_v_SH = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_SH,(/1,2/)),2),1)
   dum_v_both = dim_avg_n(dim_avg_n(dum_v,2),1)
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
   bw_data(ind_st:ind_ed) = dum_v_both
   bw_data_NH(ind_st:ind_ed) = dum_v_NH
   bw_data_SH(ind_st:ind_ed) = dum_v_SH
   bw_time(ind_st:ind_ed) = dum_t
   delete([/fin,dum_t,dum_v,dum_v_both,dum_v_NH,dum_v_SH,dim_t/])
end do

;;;;;;;;;;
; Get annual means now
;;;;;;;;;;
; Global
b_data_ann_mean = avg(b_data)
bw_data_ann_mean = avg(bw_data)
print("Global annual mean sfc temp over entire PI control runs")
print("CAM: "+b_data_ann_mean+" K")
print("WACCM: "+bw_data_ann_mean+" K")
b_data_ann_stddev = stddev(b_data)
bw_data_ann_stddev = stddev(bw_data)
print("Global annual stddev sfc temp over entire PI control runs")
print("CAM: "+b_data_ann_stddev+" K")
print("WACCM: "+bw_data_ann_stddev+" K")

; NH
b_data_NH_ann_mean = avg(b_data_NH)
bw_data_NH_ann_mean = avg(bw_data_NH)
print("NH annual mean sfc temp over entire PI control runs")
print("CAM: "+b_data_NH_ann_mean+" K")
print("WACCM: "+bw_data_NH_ann_mean+" K")
b_data_NH_ann_stddev = stddev(b_data_NH)
bw_data_NH_ann_stddev = stddev(bw_data_NH)
print("NH annual stddev sfc temp over entire PI control runs")
print("CAM: "+b_data_NH_ann_stddev+" K")
print("WACCM: "+bw_data_NH_ann_stddev+" K")

; SH
b_data_SH_ann_mean = avg(b_data_SH)
bw_data_SH_ann_mean = avg(bw_data_SH)
print("SH annual mean sfc temp over entire PI control runs")
print("CAM: "+b_data_SH_ann_mean+" K")
print("WACCM: "+bw_data_SH_ann_mean+" K")
b_data_SH_ann_stddev = stddev(b_data_SH)
bw_data_SH_ann_stddev = stddev(bw_data_SH)
print("SH annual stddev sfc temp over entire PI control runs")
print("CAM: "+b_data_SH_ann_stddev+" K")
print("WACCM: "+bw_data_SH_ann_stddev+" K")

;;;;;;;;;;;;;;;;;;;;;; END script
end
