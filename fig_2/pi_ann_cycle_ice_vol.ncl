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

b_files = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode+"/gn/latest/*.nc")
b_files_2 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode2+"/gn/latest/*.nc")
nfiles_b = dimsizes(b_files)
bw_files = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode+"/gn/latest/*.nc")
bw_files_2 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode2+"/gn/latest/*.nc")
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
kmt = b_file->KMT
; mask tarea where there is land
tarea = where(kmt .gt. 0.5, tarea, tarea@_FillValue)
; convert from cm^2 to m^2
tarea = tarea/100.0/100.0
dims = dimsizes(lat2d)
nlat = dims(0)
nlon = dims(1)
delete(b_file)

nmons = 12
nyrs_b = 1200
nyrs_bw = 499
ntimes_b = nmons*nyrs_b
ntimes_bw = nmons*nyrs_bw

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
print("Loading b expt PI control")
b_data_all = new((/ntimes_b,nlat,nlon/),"float")
b_data_2_all = new((/ntimes_b,nlat,nlon/),"float")
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   fin = addfile(b_files(n),"r")
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
   b_time(ind_st:ind_ed) = dum_t
   ; deal with variable of interest
   b_data_all(ind_st:ind_ed,:,:) = fin->$varcode$
   delete([/fin,dum_t,dim_t/])
   fin = addfile(b_files_2(n),"r")
   b_data_2_all(ind_st:ind_ed,:,:) = fin->$varcode2$
   delete(fin)
end do

;;; bw files
print("Loading bw expt PI control")
bw_data_all = new((/ntimes_bw,nlat,nlon/),"float")
bw_data_2_all = new((/ntimes_bw,nlat,nlon/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   fin = addfile(bw_files(n),"r")
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
   bw_time(ind_st:ind_ed) = dum_t
   ; deal with variable of interest
   bw_data_all(ind_st:ind_ed,:,:) = fin->$varcode$
   delete([/fin,dum_t,dim_t/])
   fin = addfile(bw_files_2(n),"r")
   bw_data_2_all(ind_st:ind_ed,:,:) = fin->$varcode2$
   delete(fin)
end do

;;;;;;;;;;
; Calculate variable over hemisphere
;;;;;;;;;;
print("calculate hemispheric volume")
; ice vol = sum of actual ice thickness per cell times the area of that cell for all grid cells with concentration > 15%
; ice vol = sum(tarea*sivol)

; calculate hemispheric sum 
b_data = dim_sum_n(dim_sum_n(b_data_all*conform(b_data_all,mask_hemi,(/1,2/)),2),1)
bw_data = dim_sum_n(dim_sum_n(bw_data_all*conform(bw_data_all,mask_hemi,(/1,2/)),2),1)

; convert to 10^3 km^3 for best y axis units
b_data = b_data/(10.^9)/(10.^3)
bw_data = bw_data/(10.^9)/(10.^3)

units = "10~S~3~N~ km~S~3~N~"

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
b_data_var = new((/nmons/),typeof(b_data))
bw_data_var = new((/nmons/),typeof(bw_data))

m = 0
do m = 0, nmons -1
  inds = ind((mm_b .eq. (m+1)).and.((yy_b .ge. 100.0).and.(yy_b .le. 499.0)))
  b_data_avg(m) = dim_avg_n(b_data(inds),0)
  b_data_stddev(m) = dim_stddev_n(b_data(inds),0)
  b_data_var(m) = dim_variance_n(b_data(inds),0)
  npts_b = dimsizes(inds)
  delete(inds)
  inds = ind((mm_bw .eq. (m+1)).and.((yy_bw .ge. 100.0).and.(yy_bw .le. 499.0)))
  bw_data_avg(m) = dim_avg_n(bw_data(inds),0)
  bw_data_stddev(m) = dim_stddev_n(bw_data(inds),0)
  bw_data_var(m) = dim_variance_n(bw_data(inds),0)
  npts_bw = dimsizes(inds)
  delete(inds)
end do

;;;;;;;;;;
; Find if means and variances are statistically different
;;;;;;;;;;
print("Checking significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
prob1 = 100.*(1. - ttest(b_data_avg,b_data_var,399,bw_data_avg,bw_data_var,399,False,False))

; use f-test to see if variances are different
prob2 = 100.*(1.-ftest(b_data_var,399,bw_data_var,399,0))

;;;;;;;;;;
; Differences for plotting
;;;;;;;;;;

data_diff = b_data_avg - bw_data_avg
stddev_diff = b_data_stddev - bw_data_stddev

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

data_sub = new((/4,nmons/),typeof(b_data))
data_sub(0,:) = b_data_avg
data_sub(1,:) = bw_data_avg
data_sub(2,:) = where(prob1 .ge. 95.0, b_data_avg, data_sub@_FillValue)
data_sub(3,:) = where(prob1 .ge. 95.0, bw_data_avg, data_sub@_FillValue)

data_stddev = new((/4,nmons/),typeof(b_data))
data_stddev(0,:) = b_data_stddev
data_stddev(1,:) = bw_data_stddev
data_stddev(2,:) = where(prob2 .ge. 95.0, b_data_stddev, data_sub@_FillValue)
data_stddev(3,:) = where(prob2 .ge. 95.0, bw_data_stddev, data_sub@_FillValue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
outdir = "./"
fout0 = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+vartitle+"_for_annual_cycle"

ncdf = addfile(fout0+".nc","c")
ncdf->data_sub = data_sub
ncdf->data_stddev = data_stddev

print("Output file with values")
return

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+vartitle)

; set outdir and fout names
   outdir = "./"

xdata = ispan(1, 12, 1)
label_inds = xdata
label_str = flt2string(int2flt(xdata))

if (hemi .eq. "NH")
   pref = "Arctic - "
   ymin = 5.0
   ymax = 38.0
   yminor = 4
   x0 = 9.0
   x1 = 3.0
   x2 = x1
   x3 = 4.0
   y0 = 34.0
   y1 = 19.0
   y2 = 17.0
   y3 = 2.15
end if
if (hemi .eq. "SH")
   pref = "Antarctic - "
   ymin = 2.0
   ymax = 24.0
   yminor = 3
   x0 = 3.0
   x1 = 3.0
   x2 = x1
   x3 = 4.0
   y0 = 22.0
   y1 = 8.0
   y2 = 7.0
   y3 = 1.02
end if

;;;;;;;;;;;;;
; XY Plot Resources
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
  ;opts@tmYLMinorPerMajor           = yminor
  ;opts@tmYRMinorPerMajor           = yminor
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@gsnYRefLine                 = 0.0

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2

;; resources for bw
  opts@tiMainOn = False
  opts@xyMarkLineMode = "MarkLines"
  opts@xyLineColors = (/"black","firebrick","black","firebrick"/)
  opts@xyLineThicknessF = 3.5
  opts@xyMarkers = (/1,1,16,16/)
  opts@xyMarkerColors = opts@xyLineColors
  opts@xyDashPatterns = (/0,0,0,0/)
  txres_b = txres
  txres_b@txFontColor = "black"
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
 fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+vartitle+"_annual_cycle"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(2, graphic)    

; seasonal cycle both
p = 0
lplot = gsn_csm_xy(wks,xdata,data_sub,opts)
dum0 = gsn_add_text(wks,lplot,"a) "+pref+"PI-Control",x0,y0,txres_b)
plot(p) = lplot
; seasonal stddev
p = 1
lplot = gsn_csm_xy(wks,xdata,data_stddev,opts)
dum3 = gsn_add_text(wks,lplot,"c) "+pref+"Standard Deviation",x3,y3,txres_b)
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

print("Figure 1 done, yay!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
