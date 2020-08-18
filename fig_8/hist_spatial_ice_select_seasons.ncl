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
;varcode = "siconc"
;hemi = "NH"
; ******************************************************
; BEGIN SCRIPT
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ssmi/"+hemi+"/"
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
delete([/b_file,dims/])

; get hemispheric mask
mask_hemi = tofloat(tarea)

if (hemi .eq. "NH")then
   mask_hemi = where(lat2d .ge. 60.0, mask_hemi, mask_hemi@_FillValue)
end if
if (hemi .eq. "SH")then
   mask_hemi = where(lat2d .le. -50.0, mask_hemi, mask_hemi@_FillValue)
end if

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
mons = fspan(1.,12.,12)
nmons = dimsizes(mons)

delete([/b_file,time_dum1,time_dum2,time_dum3,ind1,ind2/])

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data")
; data is in 50 or 100 year chunks

;;; b files
print("Loading b expt historical")
b_data = new((/nens_b,ntimes,nlat,nlon/),"float")
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
      b_data(n,:,:,:) = fin->$varcode$(ind1:ind2,:,:)
      delete([/fin,time,time_dum,time_str/])
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
      b_data(n,0:(ind2-ind1),:,:) = fin->$varcode$(ind1:ind2,:,:)
      delete([/fin,time,time_dum,time_str/])
      ; second chunk
      fin = addfile(files(3),"r")
      b_data(n,(ind2-ind1+1):ntimes-1,:,:) = fin->$varcode$(:,:,:)
      delete([/fin,ind1,ind2/])
   end if
   delete(files)
end do

;;; bw files
; loop through bw files to get the data we need
bw_data = new((/nens_bw,ntimes,nlat,nlon/),"float")
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
      bw_data(n,:,:,:) = fin->$varcode$(ind1:ind2,:,:)
      delete([/fin,time,time_dum,time_str/])
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
      bw_data(n,0:(ind2-ind1),:,:) = fin->$varcode$(ind1:ind2,:,:)
      delete([/fin,time,time_dum,time_str/])
      ; second chunk
      fin = addfile(files(3),"r")
      bw_data(n,(ind2-ind1+1):ntimes-1,:,:) = fin->$varcode$(:,:,:)
      delete([/fin,ind1,ind2/])
   end if
   delete(files)
end do

;;;;;;;;;;
; Load observed data
;;;;;;;;;;
print("Loading SSMI data")
fin = addfile(data_dir+"ssmi_monthly_data_1979-2014.nc","r")
lat2d_obs = fin->Lat
lon2d_obs = fin->Lon
dims = dimsizes(lat2d_obs)
nlat_obs = dims(0)
nlon_obs = dims(1)
delete(dims)
time_obs = fin->Time
; get subset of years/months that correspond to obs time
yy_obs = yy_b_all(48:479)
mm_obs = mm_b_all(48:479)

; Use bootstrap data based on Laura recommendation
obs_data = fin->Bootstrap
;; mask land values
obs_data@_FillValue = 9.96921e+36
obs_data = where(obs_data .ge. 200.0, obs_data@_FillValue, obs_data)
obs_data = where(obs_data .lt. 0.0, obs_data@_FillValue, obs_data)
delete(fin)

;return

;;;;;;;;;;
; Calculate means, variances, and statistical significance of difference
;;;;;;;;;;
; loop through data to get means over: ensembles, years, and months
; n = nens * nyrs * nmonths

print ("Getting seasonal averages")
nseas = 4
b_data_seas = new((/nseas,nlat,nlon/),"float")
bw_data_seas = new((/nseas,nlat,nlon/),"float")
b_data_seas_var = new((/nseas,nlat,nlon/),"float")
bw_data_seas_var = new((/nseas,nlat,nlon/),"float")
obs_data_seas = new((/nseas,nlat_obs,nlon_obs/),"float")

s = 0 
do s = 0, nseas -1
  ; b-case data
  ; JFM
  if (s .eq. 0)
     inds = ind(((mm_b_all .eq. 1) .or. (mm_b_all .eq. 2) .or. (mm_b_all .eq. 3)) .and. (yy_b_all .ge. 1979.0))
  end if
  ; AMJ
  if (s .eq. 1)
     inds = ind(((mm_b_all .eq. 4) .or. (mm_b_all .eq. 5) .or. (mm_b_all .eq. 6)) .and. (yy_b_all .ge. 1979.0))
  end if 
  ; JAS
  if (s .eq. 2)
     inds = ind(((mm_b_all .eq. 7) .or. (mm_b_all .eq. 8) .or. (mm_b_all .eq. 9)) .and. (yy_b_all .ge. 1979.0))
  end if 
  ; OND
  if (s .eq. 3)
     inds = ind(((mm_b_all .eq. 10) .or. (mm_b_all .eq. 11) .or. (mm_b_all .eq. 12)) .and. (yy_b_all .ge. 1979.0))
  end if 
     ; reshape data to make it easier to calculate over
     ntot_b = dimsizes(inds)*nens_b
     dum = reshape(b_data(:,inds,:,:),(/ntot_b,nlat,nlon/))
     b_data_seas(s,:,:) = dim_avg_n_Wrap(dum,0)
     b_data_seas_var(s,:,:) = dim_variance_n_Wrap(dum,0)
     delete([/inds,dum/])
  ; bw-case data
  ; JFM
  if (s .eq. 0)
     inds = ind(((mm_b_all .eq. 1) .or. (mm_b_all .eq. 2) .or. (mm_b_all .eq. 3)) .and. (yy_b_all .ge. 1979.0))
  end if
  ; AMJ
  if (s .eq. 1)
     inds = ind(((mm_b_all .eq. 4) .or. (mm_b_all .eq. 5) .or. (mm_b_all .eq. 6)) .and. (yy_b_all .ge. 1979.0))
  end if 
  ; JAS
  if (s .eq. 2)
     inds = ind(((mm_b_all .eq. 7) .or. (mm_b_all .eq. 8) .or. (mm_b_all .eq. 9)) .and. (yy_b_all .ge. 1979.0))
  end if 
  ; OND
  if (s .eq. 3)
     inds = ind(((mm_b_all .eq. 10) .or. (mm_b_all .eq. 11) .or. (mm_b_all .eq. 12)) .and. (yy_b_all .ge. 1979.0))
  end if
     ntot_bw = dimsizes(inds)*nens_bw
     dum = reshape(bw_data(:,inds,:,:),(/ntot_bw,nlat,nlon/))
     bw_data_seas(s,:,:) = dim_avg_n_Wrap(dum,0)
     bw_data_seas_var(s,:,:) = dim_variance_n_Wrap(dum,0)
     delete([/inds,dum/])
end do

;;;;;;;;;;
; Get statistical significance at each point
;;;;;;;;;;
print("Calculating significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level
prob = 100.*(1. - ttest(b_data_seas,b_data_seas_var,ntot_b,bw_data_seas,bw_data_seas_var,ntot_bw,True,False))

;;;;;;;;;;
; Get Averages for obs
;;;;;;;;;;
print ("Getting monthly and seasonal satellite averages")
obs_data_avg = new((/nmons,nlat_obs,nlon_obs/),"float")

m = 0
do m = 0, nmons -1
   mm = mons(m)
   ; average the obs data
   inds = ind(mm_obs .eq. mm)
   obs_data_avg(m,:,:) = dim_avg_n(obs_data(inds,:,:),0)
   ; mask observational pole hole
   obs_data_avg(m,:,:) = where(lat2d_obs .ge. 86.5, 100.0, obs_data_avg(m,:,:))
   delete(inds)
   delete(mm)
end do

; seasonal
obs_data_seas = new((/nseas,nlat_obs,nlon_obs/),"float")
;JFM average
obs_data_seas(0,:,:) = dim_avg_n(obs_data_avg(0:2,:,:),0)
;AMJ average
obs_data_seas(1,:,:) = dim_avg_n(obs_data_avg(3:5,:,:),0)
;JAS average
obs_data_seas(2,:,:) = dim_avg_n(obs_data_avg(6:8,:,:),0)
;OND average
obs_data_seas(3,:,:) = dim_avg_n(obs_data_avg(9:11,:,:),0)

;;;;;;;;;;
; Get differences
;;;;;;;;;;
; calculate difference
diff_seas = b_data_seas - bw_data_seas

; mask differences where there is no ice
diff_seas = where(b_data_seas .lt. 15.0, 0.0, diff_seas)

;;;;;;;;;;
; Assign lat/lon for plotting
;;;;;;;;;;
b_data_seas@lat2d = lat2d
b_data_seas@lon2d = lon2d
bw_data_seas@lat2d = lat2d
bw_data_seas@lon2d = lon2d
diff_seas@lat2d = lat2d
diff_seas@lon2d = lon2d
prob@lat2d = lat2d
prob@lon2d = lon2d

obs_data_avg@lat2d = lat2d_obs
obs_data_avg@lon2d = lon2d_obs
obs_data_seas@lat2d = lat2d_obs
obs_data_seas@lon2d = lon2d_obs

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set contour levels based on varcode
if (varcode .eq. "siconc")
   main_cntrs=(/15,30,40,50,60,70,80,90,95/)
   diff_cntrs=(/-25,-20,-15,-10,-5,-1,1,5,10,15,20,25/)
end if
mons_str = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct","Nov","Dec"/)
if (hemi .eq. "NH")
   seas_str = (/"Winter - JFM","Summer - JAS"/)
end if
if (hemi .eq. "SH")
   seas_str = (/"Summer - JFM","Winter - JAS"/)
end if

; set outdir and fout names
   outdir = "./"

; set general resources
  res                     = True          ; plot mods desired
  res@gsnDraw             = False         ; don't draw the plot yet
  res@gsnFrame            = False         ; or the frame
  res@cnLinesOn           = False         ; turn off contour lines
  res@txFontHeightF       = 0.030          ; set subtitle text font
  res@gsnLeftString       = ""
  res@gsnRightString      = ""
  res@gsnCenterString     = ""
  res@gsnAddCyclic        = True          ; global data, wrap around 
 
; Map resources for main dataset
  res@gsnPolar             = hemi
  if (hemi .eq. "NH")
     res@mpMinLatF            = 50.0
  end if
  if (hemi .eq. "SH")
     res@mpMaxLatF            = -50.0
  end if
  res@mpGridAndLimbOn      = False            ; Turn on lat/lon grid
  res@mpFillOn             = True             ; Enables map fill
  res@mpLandFillColor      = "grey"
  res@mpFillDrawOrder      = "PostDraw"
  res@mpPerimDrawOrder     = "PostDraw"
  res@mpGeophysicalLineColor      = "black"          ; Color of continents
  res@mpGeophysicalLineThicknessF = 1.5
  res@gsnPolarLabelSpacing = 45
  res@gsnPolarLabelFontHeightF = 0.02  

  ; contour resources for everything
  res@cnFillOn            = True          ; turn on color
  res@cnFillMode          = "AreaFill"
  res@cnLinesOn           = False
  res@lbLabelBarOn        = False         ; turn off label bars
  res@cnInfoLabelOn       = False         ; turn off contour info labels
  res@cnLineLabelsOn      = False 
  res@cnMissingValFillColor = "grey" ;white"
  res@cnLevelSelectionMode    = "ExplicitLevels"

; set resources for 0% concentration overlay
  res_0                            = True               ; reset options
  res_0@gsnLeftString              = ""
  res_0@gsnRightString             = ""
  res_0@gsnFrame                   = False             ; Do not draw plot
  res_0@gsnDraw                    = False             ; Do not advance fr
  res_0@cnInfoLabelOn              = False             ; Turn off info on contours at bottom
  res_0@lbLabelBarOn               = False             ; Turn off label bar
  res_0@cnFillOn                   = True              ; turn on filled contours
  res_0@cnMonoFillColor            = True
  res_0@cnFillColor                = "grey"
  res_0@cnMonoFillPattern          = False
  res_0@cnLevelSelectionMode       = "ManualLevels"    ; set manual contour
  res_0@cnMinLevelValF             = 0.0              ; set min contour level specified above
  res_0@cnMaxLevelValF             = 0.5              ; set max contour level specified above
  res_0@cnLevelSpacingF            = 0.5              ; set contour spacing specified above
  res_0@cnFillPatterns             = (/-1,17,-1/)      ; stipling where significant and above threshold

; set resources for 95% statistical significance overlay
  res_95                            = True               ; reset options
  res_95@gsnLeftString              = ""
  res_95@gsnRightString             = ""
  res_95@gsnFrame                   = False             ; Do not draw plot
  res_95@gsnDraw                    = False             ; Do not advance fr
  res_95@cnInfoLabelOn              = False             ; Turn off info on contours at bottom
  res_95@lbLabelBarOn               = False             ; Turn off label bar
  res_95@cnFillOn                   = True              ; turn on filled contours
  res_95@cnMonoFillColor            = True
  res_95@cnMonoFillPattern          = False
  res_95@cnLevelSelectionMode       = "ManualLevels"    ; set manual contour
  res_95@cnMinLevelValF             = 0.0              ; set min contour level specified above
  res_95@cnMaxLevelValF             = 95.              ; set max contour level specified above
  res_95@cnLevelSpacingF            = 95.0              ; set contour spacing specified above
  res_95@cnFillPatterns             = (/-1,17,-1/)      ; stipling where significant and above threshold

; resources for 15% observational overlay
  res15per = True
  res15per@gsnDraw          = False     ; don't draw
  res15per@gsnFrame         = False     ; don't advance frame
  res15per@cnFillMode       = "RasterFill"
  res15per@cnLevelSelectionMode =  "ExplicitLevels"
  res15per@cnLevels         = (/15.0/)
  res15per@cnLinesOn        = True
  res15per@cnInfoLabelOn    = False
  res15per@gsnLeftString    = ""
  res15per@gsnRightString   = ""
  res15per@tiXAxisString    = ""
  res15per@tiYAxisString    = ""
  res15per@cnLineThicknessF = 3.5 ;2.5       ; thicker contours
  res15per@cnLineColor      = "orange" ;red" ;darkgoldenrod4" 
  res15per@cnLineDashPattern = 0
  res15per@cnLineLabelsOn   = False     ; no line labels
  res15per@trGridType                 = "TriangularMesh"

; resources for line
  lnres                             = True
  lnres@gsLineThicknessF            = 0.5
  lnres@gsLineColor                 = "black"
  lnres@gsLineDashPattern           = 0

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True	       ; fill the page
  pnlres@gsnPanelLabelBar           = True
  pnlres@lbOrientation              = "vertical"

;;;;;;;;;;;;;
; Make regular plots
;;;;;;;;;;;;;
print("Plotting seasonal data")
;set info for this plot
   fout = outdir+"CMIP6_b_bw_"+hemi+"_hist_"+varcode+"_seasonal"
   title_fin = b_data@title + " ("+b_data@units+") WACCM"

; open wks
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
  ;gsn_define_colormap(wks,"julies")
  gsn_define_colormap(wks,"CBR_wet")
  gsn_reverse_colormap(wks)

; set resources
  res_main                = res
  res_main@cnLevels       = main_cntrs

; do individual plots
plot = new(2, graphic)
plot_obs = new(2,graphic)
plot_ocn = new(2,graphic)
; Winter
s = 0
res_main@gsnCenterString = seas_str(s)
cplot = gsn_csm_contour_map(wks,bw_data_seas(0,:,:),res_main)
plot_obs(s) = gsn_csm_contour(wks,obs_data_seas(0,:,:), res15per)
overlay(cplot,plot_obs(s))
plot_ocn(s) = gsn_csm_contour(wks,bw_data_seas(0,:,:),res_0)
overlay(cplot,plot_ocn(s))
plot(s) = cplot

; Summer
s = 1
res_main@gsnCenterString = seas_str(s)
cplot = gsn_csm_contour_map(wks,bw_data_seas(2,:,:),res_main)
plot_obs(s) = gsn_csm_contour(wks,obs_data_seas(2,:,:), res15per)
overlay(cplot,plot_obs(s))
plot_ocn(s) = gsn_csm_contour(wks,bw_data_seas(2,:,:),res_0)
overlay(cplot,plot_ocn(s))
plot(s) = cplot

; panel plot
print("Paneling")
  ; make panel plot
  ;pnlres@txString  = title_fin
  gsn_panel(wks,(/plot/),(/1,2/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,title_fin,plot/])

;;;;;;;;;;;;;
; Make difference plots
;;;;;;;;;;;;;
print("Plotting seasonal difference data")
;set info for this plot
   fout = outdir+"CMIP6_b_bw_"+hemi+"_hist_"+varcode+"_diff_seasonal"
   title_fin = b_data@title + " ("+b_data@units+") difference: CAM-WACCM"

; open wks
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
  ;gsn_define_colormap(wks,"GMT_polar")
  gsn_define_colormap(wks,"temp_diff_18lev")

; set resources
  res_diff                = res
  res_diff@cnLevels       = diff_cntrs

; do individual plots
plot = new(2, graphic)
plot_sig = new(2,graphic)
; Winter
s = 0
res_diff@gsnCenterString = seas_str(s)
cplot = gsn_csm_contour_map(wks,diff_seas(0,:,:),res_diff)
plot_sig(s) = gsn_csm_contour(wks,prob(0,:,:),res_95)
overlay(cplot,plot_sig(s))
plot(s) = cplot
; Summer
s = 1
res_diff@gsnCenterString = seas_str(s)
cplot = gsn_csm_contour_map(wks,diff_seas(2,:,:),res_diff)
plot_sig(s) = gsn_csm_contour(wks,prob(2,:,:),res_95)
overlay(cplot,plot_sig(s))
plot(s) = cplot

; panel plot
print("Paneling")
  ; make panel plot
  ;pnlres@txString  = title_fin
  gsn_panel(wks,(/plot/),(/1,2/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,title_fin,plot/])

print("Done with all plots. Hooray!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
