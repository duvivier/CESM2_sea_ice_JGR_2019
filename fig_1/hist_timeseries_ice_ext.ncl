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

mm = 13.0       ; Choose month 1-12 or 13, which is annual mean
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

; preallocate for files
b_files = new((/nens_b/),"string")
bw_files = new((/nens_bw/),"string")

n = 0
do n = 0, nens_b -1
   i = n+1
   root = "r"+flt2string(int2flt(i))+"i1p1f1"
   b_files(n) = cmip_dir+"/CESM2/historical/"+root+"/SImon/"+varcode+"/gn/latest/"
   delete(root)
end do

n = 0
do n = 0, nens_bw -1
   i = n+1
   root = "r"+flt2string(int2flt(i))+"i1p1f1"
   bw_files(n) = cmip_dir+"/CESM2-WACCM/historical/"+root+"/SImon/"+varcode+"/gn/latest/"
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

; get strings to be correct dates for reality
mm_b_all = mm_b_all -1
inds = ind(mm_b_all .eq. 0)
mm_b_all(inds) = 12
yy_b_all(inds) = yy_b_all(inds)-1
time_b_str = flt2string(yy_b_all)+"-"+flt2string(mm_b_all)
years_all = onedtond(yy_b_all, (/ntimes/12, 12/))
years = years_all(:,1)
nyrs = dimsizes(years)
nmons = 12

delete([/b_file,time_dum1,time_dum2,time_dum3,inds,ind1,ind2/])

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
; ice ext = grid cells with concentration > 15%, area of these cells (m^2)
; sum over just hemi cells

;;; b files
; loop through b files to get the data we need
b_data = new((/nens_b,ntimes/),"float")
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
      dum_v = fin->$varcode$(ind1:ind2,:,:)
      ; mask the var by where there is concentration > 15%,
      ; then mask just hemisphere of interest, then sum up to get single value
      b_data(n,:) = dim_sum_n(dim_sum_n(where(dum_v .ge. 15.0,1.0,0.0)*conform(dum_v,mask_hemi,(/1,2/)),2),1)  
      delete([/fin,time,time_dum,time_str,dum_v/])
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
      dum_v = fin->$varcode$(ind1:ind2,:,:)
      ; mask the var by where there is concentration > 15%,
      ; then mask just hemisphere of interest, then sum up to get single value
      b_data(n,0:(ind2-ind1)) = dim_sum_n(dim_sum_n(where(dum_v .ge. 15.0,1.0,0.0)*conform(dum_v,mask_hemi,(/1,2/)),2),1)
      delete([/fin,time,time_dum,time_str,dum_v/])
      ; second chunk
      fin = addfile(files(3),"r")
      dum_v = fin->$varcode$(:,:,:)
      b_data(n,(ind2-ind1+1):ntimes-1) = dim_sum_n(dim_sum_n(where(dum_v .ge. 15.0,1.0,0.0)*conform(dum_v,mask_hemi,(/1,2/)),2),1)
      delete([/fin,ind1,ind2,dum_v/])
   end if
   delete(files)
end do

;;; bw files
; loop through bw files to get the data we need
bw_data = new((/nens_bw,ntimes/),"float")
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
      dum_v = fin->$varcode$(ind1:ind2,:,:)
      ; mask the var by where there is concentration > 15%,
      ; then mask just hemisphere of interest, then sum up to get single value
      bw_data(n,:) = dim_sum_n(dim_sum_n(where(dum_v .ge. 15.0,1.0,0.0)*conform(dum_v,mask_hemi,(/1,2/)),2),1)
      delete([/fin,time,time_dum,time_str,dum_v/])
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
      dum_v = fin->$varcode$(ind1:ind2,:,:)
      ; mask the var by where there is concentration > 15%,
      ; then mask just hemisphere of interest, then sum up to get single value
      bw_data(n,0:(ind2-ind1)) = dim_sum_n(dim_sum_n(where(dum_v .ge. 15.0,1.0,0.0)*conform(dum_v,mask_hemi,(/1,2/)),2),1)
      delete([/fin,time,time_dum,time_str,dum_v/])
      ; second chunk
      fin = addfile(files(3),"r")
      dum_v = fin->$varcode$(:,:,:)
      bw_data(n,(ind2-ind1+1):ntimes-1) = dim_sum_n(dim_sum_n(where(dum_v .ge. 15.0,1.0,0.0)*conform(dum_v,mask_hemi,(/1,2/)),2),1)
      delete([/fin,ind1,ind2,dum_v/])
   end if
   delete(files)
end do

; convert to km2
b_data = b_data/(10.^6)
bw_data = bw_data/(10.^6)

; convert to 10^6 km^2 for better y axis units 
b_data = b_data/(10.^6)
bw_data = bw_data/(10.^6)

;;;;;;;;;;
; Get just time of interest
;;;;;;;;;;
; preallocate arrays
b_data_sub = new((/nens_b+1,nyrs/),typeof(b_data))
bw_data_sub = new((/nens_bw+1,nyrs/),typeof(bw_data))

; choose months of interest
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
   inds = ind(mm_b_all .eq. mm)
   b_data_sub(0:nens_b-1,:) = b_data(:,inds)
   bw_data_sub(0:nens_bw-1,:) = bw_data(:,inds)
   delete(inds)
end if

; for annual mean
if (mm .eq. 13.0)
   n = 0
   do n = 0, nens_b-1
      b_data_tmp = onedtond(b_data(n,:),(/nyrs,nmons/))
      b_data_sub(n,:) = dim_avg_n_Wrap(b_data_tmp,1)
      delete(b_data_tmp)
   end do

   n = 0
   do n = 0, nens_bw-1
      bw_data_tmp = onedtond(bw_data(n,:),(/nyrs,nmons/))
      bw_data_sub(n,:) = dim_avg_n_Wrap(bw_data_tmp,1)
      delete(bw_data_tmp)
   end do
end if

;;;;;;;;;;
; Check if the model experiment ensemble means are significantly different
;;;;;;;;;;
; calculate variance
b_avg = dim_avg_n_Wrap(b_data_sub(0:nens_b-1,4:39),0)
bw_avg = dim_avg_n_Wrap(bw_data_sub(0:nens_bw-1,4:39),0)
b_var = dim_variance_n_Wrap(b_data_sub(0:nens_b-1,4:39),0)
bw_var = dim_variance_n_Wrap(b_data_sub(0:nens_bw-1,4:39),0)

; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
prob = 100.*(1. - ttest(b_avg,b_var,nens_b,bw_avg,bw_var,nens_bw,True,False))

;;;;;;;;;;
; Get NSIDC extent
;;;;;;;;;;
print("Getting NSIDC data")
if (mm .lt. 13.0)
   if (hemi .eq. "NH")
      fin = asciiread(data_dir+"/NSIDC_SeaIce_extent/N_"+mm_str+"_extent_v3.0.csv",(/40,4/),"float")
   end if
   if (hemi .eq. "SH")
      fin = asciiread(data_dir+"/NSIDC_SeaIce_extent/S_"+mm_str+"_extent_v3.0.csv",(/40,4/),"float")
   end if
   years_nsidc = fin(:,0)
   nsidc_ext_hemi = fin(:,2)  ; units of 10^6 km^2
   nsidc_ext_hemi@_FillValue = -99.9
end if

if (mm .eq. 13.0)
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
   nsidc_ext_hemi = dim_avg_n_Wrap(nsidc_all,0)
end if

;;;;;;;;;;
; Array of ensemble means and obs
;;;;;;;;;;
years_sub = ispan(1979,2014,1)
all_data_sub = new((/5,dimsizes(years_sub)/),"float")

all_data_sub(0,:) = b_avg
all_data_sub(1,:) = bw_avg
all_data_sub(2,:) = nsidc_ext_hemi(0:35)
all_data_sub(3,:) = where(prob .ge. 95.0, b_avg, all_data_sub@_FillValue)
all_data_sub(4,:) = where(prob .ge. 95.0, bw_avg, all_data_sub@_FillValue)

;;;;;;;;;;
; Calculate linear trends
;;;;;;;;;;
print("Calculating trend")
rc_b = regline_stats(years_sub,all_data_sub(0,:))
eqn_b = "y = "+rc_b+" x + "+rc_b@yintercept
print("Eqn b "+eqn_b)
print("r2 = "+rc_b@r2)

rc_bw = regline_stats(years_sub,all_data_sub(1,:))
eqn_bw = "y = "+rc_bw+" x + "+rc_bw@yintercept
print("Eqn bw "+eqn_bw)
print("r2 = "+rc_bw@r2)

rc_o = regline_stats(years_sub,all_data_sub(2,:))
eqn_o = "y = "+rc_o+" x + "+rc_o@yintercept
print("Eqn NSIDC "+eqn_o)
print("r2 = "+rc_o@r2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_hist_ice_ext"

label_inds = (/1980, 1985, 1990, 1995, 2000, 2005, 2010/)
label_minor_inds = ispan(1979,2015,1)
label_str = (/"1980", "1985", "1990", "1995", "2000", "2005", "2010"/)

if ((hemi .eq. "NH") .and. (mm .eq. 13.0))
   pref = "b) " ;Arctic - "
   ymin = 8.5
   ymax = 13.5
   yminor = 1
   x0 = 1996
   x1 = 1985
   x2 = x1
   x3 = x1
   y0 = 12.75 ;13.75
   y1 = 10.5
   y2 = 10.15
   y3 = 9.8
end if
if ((hemi .eq. "SH") .and. (mm .eq. 13.0))
   pref = "b) " ;Antarctic - "
   ymin = 9.5
   ymax = 15.5 ;16.0
   yminor = 1
   x0 = 1996
   x1 = 1985
   x2 = x1
   x3 = x1
   y0 = 13.75 ;15.75
   y1 = 10.25
   y2 = 9.8
   y3 = 9.35
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
  opts@trXMinF                     = 1979
  opts@trXMaxF                     = 2014
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
  opts_o@xyLineColors = (/"black","firebrick","blue","black","firebrick"/)
  opts_o@xyLineThicknessF = 5.0
  opts_o@xyDashPattern = 0
  ;opts_o@xyMarker = 1
  opts_o@xyMarkers = (/4,4,1,16,16/)
  opts_o@xyMarkerColors = opts_o@xyLineColors
  txres_o = txres
  txres_o@txFontColor = "blue"

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
   lplot = gsn_csm_xy(wks,years,b_data_sub,opts_b)
   oplot1 = gsn_csm_xy(wks,years,bw_data_sub,opts_bw)
   overlay(lplot,oplot1)
   oplot2 = gsn_csm_xy(wks,years_sub,all_data_sub,opts_o)
   overlay(lplot,oplot2)
   ;dum0 = gsn_add_text(wks,lplot,pref+"Ice Extent Historical",x0,y0,txres_b)
   dum1 = gsn_add_text(wks,lplot,"NSIDC",x1,y1,txres_o)
   dum2 = gsn_add_text(wks,lplot,"CAM6",x2,y2,txres_b)
   dum3 = gsn_add_text(wks,lplot,"WACCM6",x3,y3,txres_bw)
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
