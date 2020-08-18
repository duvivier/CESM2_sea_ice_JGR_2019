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
hemi = "NH"     ; choose hemisphere - "NH" or "SH"
; ******************************************************
; BEGIN SCRIPT

varcode1 = "clwvi"
varcode2 = "clivi"
varcode3 = "clt"
varcode_i = "siconca"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/NSIDC_SeaIce_extent/"
grid_dir = "/glade/p/cesm/omwg/grids/"

b_files1 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/Amon/"+varcode1+"/gn/latest/*.nc")
b_files2 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/Amon/"+varcode2+"/gn/latest/*.nc")
b_files3 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/Amon/"+varcode3+"/gn/latest/*.nc")
b_files_i = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode_i+"/gn/latest/*.nc")
nfiles_b = dimsizes(b_files1)
bw_files1 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/Amon/"+varcode1+"/gn/latest/*.nc")
bw_files2 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/Amon/"+varcode2+"/gn/latest/*.nc")
bw_files3 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/Amon/"+varcode3+"/gn/latest/*.nc")
bw_files_i = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode_i+"/gn/latest/*.nc")
nfiles_bw = dimsizes(bw_files1)

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(b_files1(0),"r")
dummy2d = b_file->$varcode1$(0,:,:)
lat1d = b_file->lat
lon1d = b_file->lon
lat2d = conform(dummy2d,lat1d,0)
copy_VarMeta(dummy2d,lat2d)
lat2d@units = "degrees"
lon2d = conform(dummy2d,lon1d,1)
delete([/b_file,dummy2d/])

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
b_data_lwp = new((/ntimes_b/),"float")
b_data_iwp = new((/ntimes_b/),"float")
b_data_clt = new((/ntimes_b/),"float")
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   ; get variables of interest
   ; deal with time index
   fin1 = addfile(b_files1(n),"r")
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
   ; get ice fraction on atm grid (0-1)
   fin_i = addfile(b_files_i(n),"r")
   dum_i = fin_i->$varcode_i$
   ; water paths
   dum_1 = fin1->$varcode1$
   fin2 = addfile(b_files2(n),"r")
   dum_2 = fin2->$varcode2$
   fin3 = addfile(b_files3(n),"r")
   dum_3 = fin3->$varcode3$
   dum_v = dum_1 - dum_2
   ; get hemispheric avgs where there is sea ice
   b_data_lwp(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_hemi,(/1,2/)),2),1)
   b_data_iwp(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum_2*conform(dum_2,mask_hemi,(/1,2/)),2),1)
   b_data_clt(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum_3*conform(dum_3,mask_hemi,(/1,2/)),2),1)
   delete([/fin1,fin2,fin3,fin_i,dum_t,dum_i,dum_1,dum_2,dum_3,dum_v,dim_t/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data_lwp = new((/ntimes_bw/),"float")
bw_data_iwp = new((/ntimes_bw/),"float")
bw_data_clt = new((/ntimes_bw/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   ; get variables of interest
   ; deal with time index
   fin1 = addfile(bw_files1(n),"r")
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
   ; get ice fraction on atm grid (0-1)
   fin_i = addfile(bw_files_i(n),"r")
   dum_i = fin_i->$varcode_i$
   ; water paths
   dum_1 = fin1->$varcode1$
   fin2 = addfile(bw_files2(n),"r")
   dum_2 = fin2->$varcode2$
   fin3 = addfile(bw_files3(n),"r")
   dum_3 = fin3->$varcode3$
   dum_v = dum_1 - dum_2
   ; get hemispheric avgs where there is sea ice
   bw_data_lwp(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum_v*conform(dum_v,mask_hemi,(/1,2/)),2),1)
   bw_data_iwp(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum_2*conform(dum_2,mask_hemi,(/1,2/)),2),1)
   bw_data_clt(ind_st:ind_ed) = dim_avg_n(dim_avg_n(dum_3*conform(dum_3,mask_hemi,(/1,2/)),2),1)
   delete([/fin1,fin2,fin3,fin_i,dum_t,dum_i,dum_1,dum_2,dum_3,dum_v,dim_t/])
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
b_data_lwp_avg = new((/nmons/),typeof(b_data_lwp))
bw_data_lwp_avg = new((/nmons/),typeof(bw_data_lwp))
b_data_lwp_stddev = new((/nmons/),typeof(b_data_lwp))
bw_data_lwp_stddev = new((/nmons/),typeof(bw_data_lwp))
b_data_lwp_var = new((/nmons/),typeof(b_data_lwp))
bw_data_lwp_var = new((/nmons/),typeof(bw_data_lwp))
b_data_iwp_avg = new((/nmons/),typeof(b_data_lwp))
bw_data_iwp_avg = new((/nmons/),typeof(bw_data_lwp))
b_data_iwp_stddev = new((/nmons/),typeof(b_data_lwp))
bw_data_iwp_stddev = new((/nmons/),typeof(bw_data_lwp))
b_data_iwp_var = new((/nmons/),typeof(b_data_lwp))
bw_data_iwp_var = new((/nmons/),typeof(bw_data_lwp))
b_data_clt_avg = new((/nmons/),typeof(b_data_lwp))
bw_data_clt_avg = new((/nmons/),typeof(bw_data_lwp))
b_data_clt_stddev = new((/nmons/),typeof(b_data_lwp))
bw_data_clt_stddev = new((/nmons/),typeof(bw_data_lwp))
b_data_clt_var = new((/nmons/),typeof(b_data_lwp))
bw_data_clt_var = new((/nmons/),typeof(bw_data_lwp))

m = 0
do m = 0, nmons -1
  inds = ind((mm_b .eq. (m+1)).and.((yy_b .ge. 100.0).and.(yy_b .le. 499.0)))
  b_data_lwp_avg(m) = dim_avg_n(b_data_lwp(inds),0)
  b_data_lwp_stddev(m) = dim_stddev_n(b_data_lwp(inds),0)
  b_data_lwp_var(m) = dim_variance_n(b_data_lwp(inds),0)
  b_data_iwp_avg(m) = dim_avg_n(b_data_iwp(inds),0)
  b_data_iwp_stddev(m) = dim_stddev_n(b_data_iwp(inds),0)
  b_data_iwp_var(m) = dim_variance_n(b_data_iwp(inds),0)
  b_data_clt_avg(m) = dim_avg_n(b_data_clt(inds),0)
  b_data_clt_stddev(m) = dim_stddev_n(b_data_clt(inds),0)
  b_data_clt_var(m) = dim_variance_n(b_data_clt(inds),0)
  delete(inds)

  inds = ind((mm_bw .eq. (m+1)).and.((yy_bw .ge. 100.0).and.(yy_bw .le. 499.0)))
  bw_data_lwp_avg(m) = dim_avg_n(bw_data_lwp(inds),0)
  bw_data_lwp_stddev(m) = dim_stddev_n(bw_data_lwp(inds),0)
  bw_data_lwp_var(m) = dim_variance_n(bw_data_lwp(inds),0)
  bw_data_iwp_avg(m) = dim_avg_n(bw_data_iwp(inds),0)
  bw_data_iwp_stddev(m) = dim_stddev_n(bw_data_iwp(inds),0)
  bw_data_iwp_var(m) = dim_variance_n(bw_data_iwp(inds),0)
  bw_data_clt_avg(m) = dim_avg_n(bw_data_clt(inds),0)
  bw_data_clt_stddev(m) = dim_stddev_n(bw_data_clt(inds),0)
  bw_data_clt_var(m) = dim_variance_n(bw_data_clt(inds),0)
  delete(inds)
end do

;;;;;;;;;;
; Check if means and variances are statistically different
;;;;;;;;;;
print("Checking significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
; LWP
prob1 = 100.*(1. - ttest(b_data_lwp_avg,b_data_lwp_var,399,bw_data_lwp_avg,bw_data_lwp_var,399,True,False))
; IWP
prob2 = 100.*(1. - ttest(b_data_iwp_avg,b_data_iwp_var,399,bw_data_iwp_avg,bw_data_iwp_var,399,True,False))
; CLT
prob3 = 100.*(1. - ttest(b_data_clt_avg,b_data_clt_var,399,bw_data_clt_avg,bw_data_clt_var,399,True,False))

;print("Probablity for CLT differences being significant "+prob3)

;;;;;;;;;;
; Calculate differences
;;;;;;;;;;
print("Calculating differences")
lwp_avg_diff = b_data_lwp_avg - bw_data_lwp_avg
dum = where(bw_data_lwp_avg .eq. 0.0, bw_data_lwp_avg@_FillValue, bw_data_lwp_avg)
lwp_avg_pcnt_diff = ((b_data_lwp_avg - bw_data_lwp_avg)/dum)*100.
lwp_avg_pcnt_diff = where(ismissing(lwp_avg_pcnt_diff), 0.0, lwp_avg_pcnt_diff)
delete(dum)
iwp_avg_diff = b_data_iwp_avg - bw_data_iwp_avg
dum = where(bw_data_iwp_avg .eq. 0.0, bw_data_iwp_avg@_FillValue, bw_data_iwp_avg)
iwp_avg_pcnt_diff = ((b_data_iwp_avg - bw_data_iwp_avg)/bw_data_iwp_avg)*100.
iwp_avg_pcnt_diff = where(ismissing(iwp_avg_pcnt_diff), 0.0, iwp_avg_pcnt_diff)
delete(dum)
clt_avg_diff = b_data_clt_avg - bw_data_clt_avg
dum = where(bw_data_clt_avg .eq. 0.0, bw_data_clt_avg@_FillValue, bw_data_clt_avg)
clt_avg_pcnt_diff = ((b_data_clt_avg - bw_data_clt_avg)/dum)*100.
clt_avg_pcnt_diff = where(ismissing(clt_avg_pcnt_diff), 0.0, clt_avg_pcnt_diff)
delete(dum)

;;;;;;;;;;
; Combine for plotting
;;;;;;;;;;
diff = new((/4,nmons+1/),typeof(b_data_lwp))
diff(0,0:nmons-1) = lwp_avg_diff
diff(1,0:nmons-1) = iwp_avg_diff
diff(2,0:nmons-1) = where(prob1 .ge. 95., lwp_avg_diff, diff@_FillValue)
diff(3,0:nmons-1) = where(prob2 .ge. 95., iwp_avg_diff, diff@_FillValue)
diff(0,nmons) = lwp_avg_diff(0)
diff(1,nmons) = iwp_avg_diff(0)
diff(2,nmons) = where(prob1(0) .ge. 95., lwp_avg_diff(0), diff@_FillValue)
diff(3,nmons) = where(prob2(0) .ge. 95., iwp_avg_diff(0), diff@_FillValue)

diff_clt = new((/2,nmons+1/),typeof(clt_avg_diff))
diff_clt(0,0:nmons-1) = clt_avg_diff
diff_clt(1,0:nmons-1) = where(prob3 .ge. 95., clt_avg_diff, diff_clt@_FillValue)
diff_clt(0,nmons) = clt_avg_diff(0)
diff_clt(1,nmons) = where(prob3(0) .ge. 95., clt_avg_diff(0), diff_clt@_FillValue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting cloud var diffs")

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_cloud_diff_annual_cycle_seaice_only"

xdata = ispan(1, nmons+1, 1)
label_inds = xdata
label_str = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)

if (hemi .eq. "NH")
   pref = "b) Cloud Differences"
   ymin = -0.02
   ymax = 0.005
   yminor = 2
   x0 = 6.0
   x1 = 2.0
   x2 = x1
   x3 = 10.5
   y0 = 0.0065
   y1 = -0.013
   y2 = -0.015
   y3 = -0.013
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
  opts@trXMaxF                     = nmons+1
  opts@tiXAxisString               = "Month"
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@gsnYRefLine                 = 0.0

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; resources for budget terms
  opts_1 = opts
  opts_1@tiYAxisString               = "kg m~S~-2~N~"
  opts_1@trYMinF                     = -0.02
  opts_1@trYMaxF                     = 0.008 ;0.02
  opts_1@tiMainOn                    = False
  opts_1@xyMarkLineMode              = "MarkLines"
  opts_1@xyLineColors                = (/"blue","deepskyblue","blue","deepskyblue"/)
  opts_1@xyMarkerColors              = opts_1@xyLineColors
  opts_1@xyMarkers                   = (/4,4,16,16/)
  opts_1@xyMarkerSizeF               = 0.025
  opts_1@xyLineThicknesses           = (/4.,4.,4.,4./)
  opts_1@xyDashPattern               = 0
; resources for snow_depth
  opts_2 = opts
  opts_2@tiYAxisString               = "%"
  opts_2@trYMinF                     = -10.
  opts_2@trYMaxF                     = 4.  ;10.
  opts_2@tiMainOn                    = False
  opts_2@xyMarkLineMode              = "MarkLines"
  opts_2@xyLineColors                = (/"firebrick","firebrick"/)
  opts_2@xyMarkerColors              = opts_2@xyLineColors
  opts_2@xyMarkers                   = (/4,16/)
  opts_2@xyMarkerSizeF               = 0.025
  opts_2@xyLineThicknesses           = (/4.,4./)
  opts_2@xyDashPattern               = 0

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2
  txres@txFontColor = "black"

; text res
  txres_1 = txres
  txres_1@txFontColor = "blue"
  txres_2 = txres
  txres_2@txFontColor = "deepskyblue"
  txres_3 = txres
  txres_3@txFontColor = "firebrick"

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
   lplot = gsn_csm_xy2(wks,xdata,diff,diff_clt,opts_1,opts_2)
   dum0 = gsn_add_text(wks,lplot,pref,x0,y0,txres)
   dum1 = gsn_add_text(wks,lplot,"LWP",x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,"IWP",x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,"CLD FRAC",x3,y3,txres_3)
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

print("Figure 1 done, yay!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
