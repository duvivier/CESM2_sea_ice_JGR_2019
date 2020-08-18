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
varcode = "net_mass_budget"
hemi = "NH"     ; choose hemisphere - "NH" or "SH"
; ******************************************************
; BEGIN SCRIPT

; load in all the variables
varcode1 = "sndmasssnf"
varcode2 = "sndmassmelt"
varcode3 = "sndmasssubl"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/NSIDC_SeaIce_extent/"
grid_dir = "/glade/p/cesm/omwg/grids/"

; original
b_files1 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode1+"/gn/latest/*.nc")
b_files2 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode2+"/gn/latest/*.nc")
b_files3 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode3+"/gn/latest/*.nc")
nfiles_b = dimsizes(b_files1)
bw_files1 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode1+"/gn/latest/*.nc")
bw_files2 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode2+"/gn/latest/*.nc")
bw_files3 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode3+"/gn/latest/*.nc")
nfiles_bw = dimsizes(bw_files1)

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

; get hemispheric mask
mask_hemi = tofloat(lat2d)

if (hemi .eq. "NH")then
   mask_hemi = where(lat2d .ge. 70.0, 1.0, mask_hemi@_FillValue)
end if
if (hemi .eq. "SH")then
   mask_hemi = where(lat2d .le. -50.0, 1.0, mask_hemi@_FillValue)
end if

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
; avg over just hemi cells

;;; b files
print("Loading b expt PI control")
b_data_net = new((/ntimes_b/),"float")
b_data_1 = new((/ntimes_b/),"float")
b_data_2 = new((/ntimes_b/),"float")
b_data_3 = new((/ntimes_b/),"float")
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   fin1 = addfile(b_files1(n),"r")
   ; deal with time index
   dum_t = fin1->time
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
   ; individual variables
   ; fin1 loaded above
   dum1 = (fin1->$varcode1$)
   b_data_1(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum1*conform(dum1,mask_hemi,(/1,2/)),2),1)
   fin2 = addfile(b_files2(n),"r")
   dum2 = (fin2->$varcode2$)
   b_data_2(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum2*conform(dum2,mask_hemi,(/1,2/)),2),1)
   fin3 = addfile(b_files3(n),"r")
   dum3 = (fin3->$varcode3$)
   b_data_3(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum3*conform(dum3,mask_hemi,(/1,2/)),2),1)
   ; get net snow mass budget
   dum_v = dum1 - dum2 - dum3 
   dum_v_hemi = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_hemi,(/1,2/)),2),1) 
   b_data_net(ind_st:ind_ed) = dum_v_hemi
   delete([/fin1,dum1,fin2,dum2,fin3,dum3,dum_t,dum_v,dum_v_hemi/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data_net = new((/ntimes_bw/),"float")
bw_data_1 = new((/ntimes_bw/),"float")
bw_data_2 = new((/ntimes_bw/),"float")
bw_data_3 = new((/ntimes_bw/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   fin1 = addfile(bw_files1(n),"r")
   ; deal with time index
   dum_t = fin1->time
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
   ; individual variables
   ; fin1 loaded above
   dum1 = (fin1->$varcode1$)
   bw_data_1(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum1*conform(dum1,mask_hemi,(/1,2/)),2),1)
   fin2 = addfile(bw_files2(n),"r")
   dum2 = (fin2->$varcode2$)
   bw_data_2(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum2*conform(dum2,mask_hemi,(/1,2/)),2),1)
   fin3 = addfile(bw_files3(n),"r")
   dum3 = (fin3->$varcode3$)
   bw_data_3(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum3*conform(dum3,mask_hemi,(/1,2/)),2),1)
   ; get net snow mass budget
   dum_v = dum1 - dum2 - dum3
   dum_v_hemi = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_hemi,(/1,2/)),2),1) 
   bw_data_net(ind_st:ind_ed) = dum_v_hemi
   delete([/fin1,dum1,fin2,dum2,fin3,dum3,dum_t,dum_v,dum_v_hemi/])
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
; Get climatologies
;;;;;;;;;;
print ("Getting annual climatology")
b_data_avg = new((/nmons/),typeof(b_data_1))
b_data_1_avg = new((/nmons/),typeof(b_data_1))
b_data_2_avg = new((/nmons/),typeof(b_data_1))
b_data_3_avg = new((/nmons/),typeof(b_data_1))
bw_data_avg = new((/nmons/),typeof(bw_data_1))
bw_data_1_avg = new((/nmons/),typeof(bw_data_1))
bw_data_2_avg = new((/nmons/),typeof(bw_data_1))
bw_data_3_avg = new((/nmons/),typeof(bw_data_1))
; variances
b_data_var = new((/nmons/),typeof(b_data_1))
b_data_1_var = new((/nmons/),typeof(b_data_1))
b_data_2_var = new((/nmons/),typeof(b_data_1))
b_data_3_var = new((/nmons/),typeof(b_data_1))
bw_data_var = new((/nmons/),typeof(bw_data_1))
bw_data_1_var = new((/nmons/),typeof(bw_data_1))
bw_data_2_var = new((/nmons/),typeof(bw_data_1))
bw_data_3_var = new((/nmons/),typeof(bw_data_1))

m = 0
do m = 0, nmons -1
  inds = ind((mm_b .eq. (m+1)).and.((yy_b .ge. 100.0).and.(yy_b .le. 499.0)))
  b_data_avg(m) = dim_avg_n(b_data_net(inds),0)
  b_data_1_avg(m) = dim_avg_n(b_data_1(inds),0)
  b_data_2_avg(m) = dim_avg_n(b_data_2(inds),0)
  b_data_3_avg(m) = dim_avg_n(b_data_3(inds),0)
  b_data_var(m) = dim_variance_n(b_data_net(inds),0)
  b_data_1_var(m) = dim_variance_n(b_data_1(inds),0)
  b_data_2_var(m) = dim_variance_n(b_data_2(inds),0)
  b_data_3_var(m) = dim_variance_n(b_data_3(inds),0)
  delete(inds)
  inds = ind((mm_bw .eq. (m+1)).and.((yy_bw .ge. 100.0).and.(yy_bw .le. 499.0)))
  bw_data_avg(m) = dim_avg_n(bw_data_net(inds),0)
  bw_data_1_avg(m) = dim_avg_n(bw_data_1(inds),0)
  bw_data_2_avg(m) = dim_avg_n(bw_data_2(inds),0)
  bw_data_3_avg(m) = dim_avg_n(bw_data_3(inds),0)
  bw_data_var(m) = dim_variance_n(bw_data_net(inds),0)
  bw_data_1_var(m) = dim_variance_n(bw_data_1(inds),0)
  bw_data_2_var(m) = dim_variance_n(bw_data_2(inds),0)
  bw_data_3_var(m) = dim_variance_n(bw_data_3(inds),0)
  delete(inds)
end do

;;;;;;;;;;
; Calculate statistical significance
;;;;;;;;;;
print("Calculating significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
prob0 = 100.*(1. - ttest(b_data_avg,b_data_var,399,bw_data_avg,bw_data_var,399,True,False))
prob1 = 100.*(1. - ttest(b_data_1_avg,b_data_1_var,399,bw_data_1_avg,bw_data_1_var,399,True,False))
prob2 = 100.*(1. - ttest(b_data_2_avg,b_data_2_var,399,bw_data_2_avg,bw_data_2_var,399,True,False))
prob3 = 100.*(1. - ttest(b_data_3_avg,b_data_3_var,399,bw_data_3_avg,bw_data_3_var,399,True,False))

;;;;;;;;;;
; Combine for plotting
;;;;;;;;;;
b_data_all = new((/4,nmons/),typeof(b_data_1))
bw_data_all = new((/4,nmons/),typeof(bw_data_1))

b_data_all(0,:) = b_data_avg
b_data_all(1,:) = b_data_1_avg
b_data_all(2,:) = -1.0*b_data_2_avg
b_data_all(3,:) = -1.0*b_data_3_avg

bw_data_all(0,:) = bw_data_avg
bw_data_all(1,:) = bw_data_1_avg
bw_data_all(2,:) = -1.0*bw_data_2_avg
bw_data_all(3,:) = -1.0*bw_data_3_avg

;;;;;;;;;;
; convert all to cm/day
;;;;;;;;;;
; snow density - 333 kg/m3
rhos = 333.0
b_data_all = (b_data_all/rhos)*100.0*60.0*60.0*24.0
bw_data_all = (bw_data_all/rhos)*100.0*60.0*60.0*24.0

;;;;;;;;;;
; Get differences
;;;;;;;;;;
diff = b_data_all - bw_data_all

diff_plots = new((/8,nmons/),typeof(b_data_1))
diff_plots(0:3,:) = diff(0:3,:)
diff_plots(4,:) = where(prob0 .ge. 95., diff_plots(0,:), diff_plots@_FillValue)
diff_plots(5,:) = where(prob1 .ge. 95., diff_plots(1,:), diff_plots@_FillValue)
diff_plots(6,:) = where(prob2 .ge. 95., diff_plots(2,:), diff_plots@_FillValue)
diff_plots(7,:) = where(prob3 .ge. 95., diff_plots(3,:), diff_plots@_FillValue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together with other data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
outdir = "./"
fout0 = outdir+"CMIP6_b_bw_"+hemi+"_piControl_snow_mass_budget_terms_for_plot"

ncdf = addfile(fout0+".nc","c")
ncdf->b_data_all          = b_data_all
ncdf->bw_data_all         = bw_data_all
ncdf->diff_plots          = diff_plots

return
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_snow_mass_budget_diffs"

xdata = ispan(1, 12, 1)
label_inds = xdata
label_str = flt2string(int2flt(xdata))

if (hemi .eq. "NH")
   pref = "b) Arctic Mass Budget"
   ymin = -0.18
   ymax = 0.14
   yminor = 4
   x0 = 3.5
   x1 = 6.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 9.0 
   x6 = x5
   x7 = x5
   x8 = x5
   x9 = x5
   y0 = 0.11
   y1 = 0.1
   y2 = 0.08
   y3 = 0.06
   y4 = 0.04
   y5 = 0.12
   y6 = 0.1
   y7 = 0.08
   y8 = 0.04
   y9 = 0.06
end if
if (hemi .eq. "SH")
   pref = "c) Antarctic Mass Budget"
   ymin = -0.05
   ymax = 0.05
   yminor = 3
   x0 = 3.5
   x1 = 7.0
   x2 = x1
   x3 = x1
   x4 = x1
   x5 = 10.0
   x6 = x5
   x7 = x5
   x8 = x5
   x9 = x5
   y0 = 0.04
   y1 = 0.045
   y2 = 0.039
   y3 = 0.033
   y4 = 0.027
   y5 = 0.045
   y6 = 0.039
   y7 = 0.033
   y8 = 0.021
   y9 = 0.027
end if


units = "cm/day"
;units = "kg m~S~-2~N~ s~S~-1~N~"

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
  opts@gsnYRefLine                 = 0.0

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; text resources
  txres = True
  txres@txFontHeightF = 0.03 ;0.02
  txres@txFontThicknessF = 2

; resources for b
  opts_b = opts
  opts_b@tiMainOn = False
  opts_b@xyMarkLineMode = "MarkLines"
  opts_b@xyLineColors = (/"black","blue","deepskyblue","purple","firebrick","palevioletred1","darkorange","green3","gold","black","blue","deepskyblue","purple","firebrick","palevioletred1","darkorange","green3","gold"/)
  opts_b@xyMarkerColors = opts_b@xyLineColors
  opts_b@xyMarkers = (/1,1,1,1,1,1,1,1,1,16,16,16,16,16,16,16,16,16/)
  opts_b@xyLineThicknesses = (/5.,4.,4.,4.,4.,4.,4.,4.,4.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
  opts_b@xyDashPattern = 0

; text
  title1 = "Net"
  txres_1 = txres
  txres_1@txFontColor = "black"
  title2 = "Frazil Growth"
  txres_2 = txres
  txres_2@txFontColor = "blue"
  title3 = "Congelation Growth"
  txres_3 = txres
  txres_3@txFontColor = "deepskyblue"
  title4 = "Snowice Growth"
  txres_4 = txres
  txres_4@txFontColor = "purple"
  title5 = "Top Melt"
  txres_5 = txres
  txres_5@txFontColor = "firebrick"
  title6 = "Bottom Melt"
  txres_6 = txres
  txres_6@txFontColor = "palevioletred1"
  title7 = "Lateral Melt"
  txres_7 = txres
  txres_7@txFontColor = "darkorange"
  title8 = "Dynamics"
  txres_8 = txres
  txres_8@txFontColor = "green3"
  title9 = "Evap/Subl"
  txres_9 = txres
  txres_9@txFontColor = "gold"

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
   lplot = gsn_csm_xy(wks,xdata,diff_plots,opts_b)
   dum0 = gsn_add_text(wks,lplot,pref+"~C~(CAM6-WACCM6)",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
   dum5 = gsn_add_text(wks,lplot,title5,x5,y5,txres_5)
   dum6 = gsn_add_text(wks,lplot,title6,x6,y6,txres_6)
   dum7 = gsn_add_text(wks,lplot,title7,x7,y7,txres_7)
   dum8 = gsn_add_text(wks,lplot,title8,x8,y8,txres_8)
   dum9 = gsn_add_text(wks,lplot,title9,x9,y9,txres_9)
   plot(0) = lplot

;;;;;;;;;;;;;
; Final Plots
;;;;;;;;;;;;;
  ; make panel plot
  ;pnlres@txString  = hemi+" Monthly Sea Ice Extent (mean over entire PIctrl)"
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Figure done, yay!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
