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
;hemi = "SH"     ; choose hemisphere - "NH" or "SH"

mm = 13.0       ; Choose month 1-12 or 13, which is annual mean
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/"
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
; ice vol = hi*tarea  ; hi is thickness per unit area, so already accounts for aice
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
   dum_v_hemi = dim_sum_n(dim_sum_n(dum_v*conform(dum_v,mask_hemi,(/1,2/)),2),1) 
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
   dum_v_hemi = dim_sum_n(dim_sum_n(dum_v*conform(dum_v,mask_hemi,(/1,2/)),2),1) 
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

; convert to km3
b_data = b_data/(10.^9)
bw_data = bw_data/(10.^9)

; convert to 10^3 km^3 for better y axis units 
b_data = b_data/(10.^3)
bw_data = bw_data/(10.^3)

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

;;;;;;;;;;
; Get running means
;;;;;;;;;;
nyrsavg = 10
b_data_ave = runave(b_data_sub,nyrsavg,0)
bw_data_ave = runave(bw_data_sub,nyrsavg,0)

;;;;;;;;;;
; Get total means and standard deviations
;;;;;;;;;;
; since we're omitting the first 100 years, take into account
b_mean = dim_avg_n(b_data_sub(99:1199),0)
b_stddev = dim_stddev_n(b_data_sub(99:1199),0)
b_var = dim_variance_n(b_data_sub(99:1199),0)

bw_mean = dim_avg_n(bw_data_sub(99:498),0)
bw_stddev = dim_stddev_n(bw_data_sub(99:498),0)
bw_var = dim_variance_n(bw_data_sub(99:498),0)

; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
prob1 = 100.*(1. - ttest(b_mean,b_var,1100,bw_mean,bw_var,399,False,False))
if (prob1 .ge. 95.0)
  sig_value1 = "SIGNIFICANT!"
else
  sig_value1 = "not significant"
end if

; use f-test to see if variances are different
prob2 = 100.*(1.-ftest(b_var,1100,bw_var,399,0))
if (prob2 .ge. 95.0)
  sig_value2 = "SIGNIFICANT!"
else
  sig_value2 = "not significant"
end if

print("Stats for entire records")
print("b_mean: "+b_mean)
print("b_stddev: "+b_stddev)
print("bw_mean: "+bw_mean)
print("bw_stddev: "+bw_stddev)
print("mean difference is "+sig_value1)
print("difference in variance is "+sig_value2)

; Now just look at overlapping times
b_mean_a = dim_avg_n(b_data_sub(99:498),0)
b_stddev_a = dim_stddev_n(b_data_sub(99:498),0)
b_var_a = dim_variance_n(b_data_sub(99:498),0)

bw_mean_a = dim_avg_n(bw_data_sub(99:498),0)
bw_stddev_a = dim_stddev_n(bw_data_sub(99:498),0)
bw_var_a = dim_variance_n(bw_data_sub(99:498),0)

; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
prob1_a = 100.*(1. - ttest(b_mean_a,b_var_a,399,bw_mean_a,bw_var_a,399,True,False))
if (prob1_a .ge. 95.0)
  sig_value1 = "SIGNIFICANT!"
else
  sig_value1 = "not significant"
end if

; use f-test to see if variances are different
prob2_a = 100.*(1.-ftest(b_var_a,399,bw_var_a,399,0))
if (prob2_a .ge. 95.0)
  sig_value2 = "SIGNIFICANT!"
else
  sig_value2 = "not significant"
end if

print("Stats for just overlapping records")
print("b_mean: "+b_mean_a)
print("b_stddev: "+b_stddev_a)
print("bw_mean: "+bw_mean_a)
print("bw_stddev: "+bw_stddev_a)
print("mean difference is "+sig_value1)
print("difference in variance is "+sig_value2)

;;;;;;;;;;
; Combine to one array
;;;;;;;;;; 
b_data_all = new((/2,nyrs_b/),"float")
bw_data_all = new((/2,nyrs_bw/),"float")

b_data_all(0,:) = b_data_sub
b_data_all(1,:) = b_data_ave

bw_data_all(0,:) = bw_data_sub
bw_data_all(1,:) = bw_data_ave

;;;;;;;;;;
; Calculate linear trend
;;;;;;;;;; 
print("Calculating trend")
rc_b = regline_stats(yy_b_sub(99:1199),b_data_sub(99:1199))
eqn_b = "y = "+rc_b+" x + "+rc_b@yintercept
print("Eqn b "+eqn_b)
print("r2 = "+rc_b@r2)
 
rc_bw = regline_stats(yy_bw_sub(99:498),bw_data_sub(99:498))
eqn_bw = "y = "+rc_bw+" x + "+rc_bw@yintercept
print("Eqn bw "+eqn_bw)
print("r2 = "+rc_bw@r2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_ice_vol_smoothed"

label_inds = (yy_b_sub(1:nyrs_b-1:100))-2
label_minor_inds = label_inds
label_str = flt2string(label_inds)

if ((hemi .eq. "NH") .and. (mm .eq. 13.0))
   pref = "c) " ;Arctic - "
   ymin = 8.0
   ymax = 32.0
   yminor = 3
   x0 = 600
   x1 = 800
   x2 = x1
   y0 = 31.0
   y1 = 30.0
   y2 = 29.0
end if
if ((hemi .eq. "SH") .and. (mm .eq. 13.0))
   pref = "c) " ;Antarctic - "
   ymin = 8.0
   ymax = 18.0
   yminor = 1
   x0 = 600
   x1 = 800
   x2 = x1
   y0 = 17.5
   y1 = 18.0
   y2 = 17.5
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
  opts@tiYAxisString               = "10~S~3~N~ km~S~3~N~"
  opts@trYMinF                     = ymin
  opts@trYMaxF                     = ymax
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@tmYLMinorPerMajor           = yminor
  opts@tmYRMinorPerMajor           = yminor
  ;opts@gsnYRefLine                 = (/b_mean,b_mean+b_stddev,b_mean-b_stddev,bw_mean,bw_mean+bw_stddev,bw_mean-bw_stddev/)
  ;opts@gsnYRefLineColors           = (/"black","black","black","firebrick","firebrick","firebrick"/)
  ;opts@gsnYRefLineThicknesses      = (/0.5,0.5,0.5,0.5,0.5,0.5/)
  ;opts@gsnYRefLineDashPatterns     = (/0,2,2,0,2,2/)

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2

; marker resources
  mkres = True
  mkres@gsMarkerSizeF  = 0.015
; marker res for means
  mkres1a = mkres
  mkres1a@gsMarkerColor = "black"
  mkres1a@gsMarkerIndex = 16
  mkres1b = mkres
  mkres1b@gsMarkerColor = "firebrick"
  mkres1b@gsMarkerIndex = 16
; marker res for stddev
  mkres2a = mkres
  mkres2a@gsMarkerIndex = 5
  mkres2a@gsMarkerThicknessF = 4.0
  mkres2a@gsMarkerColor = "black"
  mkres2b = mkres
  mkres2b@gsMarkerIndex = 5
  mkres2b@gsMarkerThicknessF = 4.0
  mkres2b@gsMarkerColor = "firebrick"

; resources for b
  opts_b = opts
  opts_b@tiMainOn = False
  opts_b@xyMarkLineMode = "Lines"
  opts_b@xyLineColors = (/"grey","black"/)
  opts_b@xyLineThicknessF = 3.5
  opts_b@xyDashPattern = 0
  txres_b = txres
  txres_b@txFontColor = "black"

; resources for bw
  opts_bw = opts
  opts_bw@tiMainOn = False
  opts_bw@xyMarkLineMode = "Lines"
  opts_bw@xyLineColors = (/"pink","firebrick"/)
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
   lplot = gsn_csm_xy(wks,yy_b_sub,b_data_all,opts_b)
   oplot1 = gsn_csm_xy(wks,yy_bw_sub,bw_data_all,opts_bw)
   overlay(lplot,oplot1)
   ;dum1 = gsn_add_text(wks,lplot,pref+"Ice Volume PI",x0,y0,txres_b)
   ;dum2 = gsn_add_text(wks,lplot,"CAM",x1,y1,txres_b)
   ;dum3 = gsn_add_text(wks,lplot,"WACCM",x2,y2,txres_bw)
   ;dum4 = gsn_add_polymarker(wks,lplot,101.0,b_mean,mkres1a)
   ;dum5 = gsn_add_polymarker(wks,lplot,1199.0,b_mean,mkres1a)
   ;dum6 = gsn_add_polymarker(wks,lplot,101.0,b_mean+b_stddev,mkres2a)
   ;dum7 = gsn_add_polymarker(wks,lplot,1199.0,b_mean+b_stddev,mkres2a)
   ;dum8 = gsn_add_polymarker(wks,lplot,101.0,b_mean-b_stddev,mkres2a)
   ;dum9 = gsn_add_polymarker(wks,lplot,1199.0,b_mean-b_stddev,mkres2a)
   ;dum10 = gsn_add_polymarker(wks,lplot,101.0,bw_mean,mkres1b)
   ;dum11 = gsn_add_polymarker(wks,lplot,1199.0,bw_mean,mkres1b)
   ;dum12 = gsn_add_polymarker(wks,lplot,101.0,bw_mean+bw_stddev,mkres2b)
   ;dum13 = gsn_add_polymarker(wks,lplot,1199.0,bw_mean+bw_stddev,mkres2b)
   ;dum14 = gsn_add_polymarker(wks,lplot,101.0,bw_mean-bw_stddev,mkres2b)
   ;dum15 = gsn_add_polymarker(wks,lplot,1199.0,bw_mean-bw_stddev,mkres2b)
   plot(0) = lplot

;;;;;;;;;;;;;
; Final Plots
;;;;;;;;;;;;;
  ; make panel plot
  ;pnlres@txString  = hemi_str+" Annual Mean Ice Volume Over Historical Period"
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Figure done, yay!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
