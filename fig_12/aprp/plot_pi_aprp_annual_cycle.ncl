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
hemi = "NH"
styr = "0400"
; ******************************************************
; BEGIN SCRIPT
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "./"
grid_dir = "/glade/p/cesm/pcwg/timeseries-cmip6/"

; case 1 - needed for grid information only
case1 = "b.e21.B1850.f09_g17.CMIP6-piControl.001"
b_files = systemfunc("ls "+grid_dir+case1+"/atm/proc/tseries/month_1/*.CLDHGH.0100*.nc")

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(b_files(0),"r")
dummy2d = b_file->CLDHGH(0,:,:)
lat1d = b_file->lat
lon1d = b_file->lon
nlat = dimsizes(lat1d)
nlon = dimsizes(lon1d)
time = b_file->time
time_conv = cd_calendar(time,0)
mm_cesm = time_conv(:,1)  ; these are off by one because files written at the end of the month!
mm_cesm = mm_cesm-1
mm_cesm = where(mm_cesm .eq. 0.0, 12.0, mm_cesm)
; now just grab part of the array (corresponds to # timesteps from APRP analysis)
mm_cesm_sub = mm_cesm(0:539)
mons = ispan(1, 12, 1)
nmons = dimsizes(mons)
delete([/b_file,dummy2d/])

;;;;;;;;;;
; Load APRP monthly Data
;;;;;;;;;;
print("Loading CESM APRP data")
; and deal with nans
dum = new((/1/),float)
fill = dum@_FillValue

fin = addfile(data_dir+"APRP_CMIP6_styr_"+styr+"_diffs_monthly.nc","r")
SW_tot = fin->monthly_SWfeedback_tot
SW_tot@_FillValue = fill
SW_tot = where(isnan_ieee(SW_tot), SW_tot@_FillValue, SW_tot)
SW_surface = fin->monthly_SWfeedback_surface
SW_surface@_FillValue = fill
SW_surface = where(isnan_ieee(SW_surface), SW_surface@_FillValue, SW_surface)
SW_cloud = fin->monthly_SWfeedback_cloud
SW_cloud@_FillValue = fill
SW_cloud = where(isnan_ieee(SW_cloud), SW_cloud@_FillValue, SW_cloud)
SW_noncloud = fin->monthly_SWfeedback_noncloud
SW_noncloud@_FillValue = fill
SW_noncloud = where(isnan_ieee(SW_noncloud), SW_noncloud@_FillValue, SW_noncloud)
delete(fin)

;;;;;;;;;;
; Get monthly means
;;;;;;;;;;
print ("Getting monthly averages")
SW_tot_avg = new((/nmons,nlat,nlon/),"float")
SW_surface_avg = new((/nmons,nlat,nlon/),"float")
SW_cloud_avg = new((/nmons,nlat,nlon/),"float")
SW_noncloud_avg = new((/nmons,nlat,nlon/),"float")
SW_tot_7090 = new((/nmons/),"float")
SW_surface_7090 = new((/nmons/),"float")
SW_cloud_7090 = new((/nmons/),"float")
SW_noncloud_7090 = new((/nmons/),"float")

inds_lat = ind(lat1d .ge. 70.0)

m = 0
do m = 0, nmons -1
   ; get data
   inds = ind(mm_cesm_sub .eq. (m+1))
   ; get monthly averages at all points
   SW_tot_avg(m,:,:) = dim_avg_n(SW_tot(inds,:,:),0)
   SW_surface_avg(m,:,:) = dim_avg_n(SW_surface(inds,:,:),0)
   SW_cloud_avg(m,:,:) = dim_avg_n(SW_cloud(inds,:,:),0)
   SW_noncloud_avg(m,:,:) = dim_avg_n(SW_noncloud(inds,:,:),0)
   ; get monthly averages just north of 70-90N
   SW_tot_7090(m) = avg(SW_tot_avg(m,inds_lat,:))
   SW_surface_7090(m) = avg(SW_surface_avg(m,inds_lat,:))
   SW_cloud_7090(m) = avg(SW_cloud_avg(m,inds_lat,:))
   SW_noncloud_7090(m) = avg(SW_noncloud_avg(m,inds_lat,:))
   delete(inds)
end do

;;;;;;;;;;;
;; Assign lat/lon for plotting
;;;;;;;;;;;
data_all = new((/4,nmons+1/),float)
data_all(0,0:nmons-1) = SW_tot_7090
data_all(1,0:nmons-1) = SW_surface_7090
data_all(2,0:nmons-1) = SW_cloud_7090
data_all(3,0:nmons-1) = SW_noncloud_7090
data_all(0,nmons) = SW_tot_7090(0)
data_all(1,nmons) = SW_surface_7090(0)
data_all(2,nmons) = SW_cloud_7090(0)
data_all(3,nmons) = SW_noncloud_7090(0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting APRP annual cycle")

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_aprp_diff_annual_cycle"

xdata = ispan(1, nmons+1, 1)
label_inds = xdata
label_str = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)

if (hemi .eq. "NH")
   pref = "a) APRP Differences"
   ymin = -3.0
   ymax = 6.0
   yminor = 3
   x0 = 3.0
   x1 = x0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 5.5
   y1 = 4.0
   y2 = 3.5
   y3 = 3.0
   y4 = 2.5
end if

units = "W m~S~-2~N~"

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
  opts@tiYAxisString               = units
  opts@trYMinF                     = ymin
  opts@trYMaxF                     = ymax
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  opts@tmYLMinorPerMajor           = yminor
  opts@tmYRMinorPerMajor           = yminor
  opts@gsnYRefLine                 = 0.0

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2
  txres@txFontColor = "black"

; resources for b
  opts = opts
  opts@tiMainOn = False
  opts@xyMarkLineMode = "MarkLines"
  opts@xyLineColors = (/"black","mediumpurple2","blue", "deepskyblue"/)
  opts@xyLineThicknessF = 5.0 ;3.5
  opts@xyDashPattern = 0
  opts@xyMarkerColors = opts@xyLineColors
  opts@xyMarkers = (/1,1,1,1/)
  txres_1 = txres
  txres_1@txFontColor = "black"
  txres_2 = txres
  txres_2@txFontColor = "mediumpurple2"
  txres_3 = txres
  txres_3@txFontColor = "blue"
  txres_4 = txres
  txres_4@txFontColor = "deepskyblue"

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
   lplot = gsn_csm_xy(wks,xdata,data_all,opts)
   dum0 = gsn_add_text(wks,lplot,pref,x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,"Total",x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,"Surface Term",x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,"Cloud Term",x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,"Noncloud Term",x4,y4,txres_4)
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

print("Done with all plots. Hooray!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
