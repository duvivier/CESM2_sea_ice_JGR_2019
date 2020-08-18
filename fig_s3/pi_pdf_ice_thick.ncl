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
;varcode = "sithick" 
;hemi = "NH"
varcode2 = "siconc"
; ******************************************************
; BEGIN SCRIPT
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/data/CMIP6/CMIP/NCAR/"
data_dir = "/glade/p/cesm/pcwg/ice/data/"
grid_dir = "/glade/p/cesm/omwg/grids/"
; set outdir and fout names
   outdir = "./"

b_files = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode+"/gn/latest/*.nc")
b_files_2 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode2+"/gn/latest/*.nc")
nfiles_b = dimsizes(b_files)
bw_files = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode+"/gn/latest/*.nc")
bw_files_2 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode2+"/gn/latest/*.nc")
nfiles_bw = dimsizes(bw_files)

;b_files = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode+"/gn/latest/*01000*.nc")
;b_files_2 = systemfunc("ls "+cmip_dir+"/CESM2/piControl/r*/SImon/"+varcode2+"/gn/latest/*01000*.nc")
;nfiles_b = dimsizes(b_files)
;bw_files = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode+"/gn/latest/*01000*.nc")
;bw_files_2 = systemfunc("ls "+cmip_dir+"/CESM2-WACCM/piControl/r*/SImon/"+varcode2+"/gn/latest/*01000*.nc")
;nfiles_bw = dimsizes(bw_files)

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

mons = fspan(1.,12.,12)
nmons = dimsizes(mons)
nyrs_b = 1200
nyrs_bw = 499
ntimes_b = nmons*nyrs_b
ntimes_bw = nmons*nyrs_bw

;;;;;;;;;;
; Get CESM mask for where there is IceSat data
;;;;;;;;;;
; mask for hemisphere
if (hemi .eq. "NH")
  ; 3 ways to choose from
  ; 1) mask of just northern pole
  ;mask_hemi = tofloat(tarea)
  ;mask_hemi = where(lat2d .ge. 70.0, 1.0, mask_hemi@_FillValue)
  ; 2) mask of whole hemi - doesn't work well
  ;fin = addfile(data_dir+"REGION_MASK_gx1v6.nc","r")
  ;mask_hemi = fin->nh_mask(0,:,:)
  ;delete(fin)
  ; 3) Mask for IceSat region only
  fin = addfile(data_dir+"CCSM4_Arctic_mask.nc","r")
  mask_hemi = fin->CCSM4_Arctic_mask_no_Barents_Kara_CAA_Nares
  delete(fin)
end if
if (hemi .eq. "SH")
  ; 2 ways to choose from
  ; 1) mask of just southern polar region
  mask_hemi = tofloat(tarea)
  mask_hemi = where(lat2d .le. -50.0, 1.0, mask_hemi@_FillValue)
  ; 2) mask of whole hemi - equal to above
  ;fin = addfile(data_dir+"REGION_MASK_gx1v6.nc","r")
  ;mask_hemi = fin->sh_mask(0,:,:)
  ;delete(fin)
end if

tarea_mask_hemi = tarea*mask_hemi

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data")
; data is in 50 or 100 year chunks

;;; b files
print("Loading b expt PI control")
b_data = new((/ntimes_b,nlat,nlon/),"float")
b_data_2 = new((/ntimes_b,nlat,nlon/),"float")
b_time = new((/ntimes_b/),"double")

n = 0
do n = 0, nfiles_b -1
   fin = addfile(b_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
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
   b_data(ind_st:ind_ed,:,:) = dum_v
   b_time(ind_st:ind_ed) = dum_t
   delete([/fin,dum_t,dum_v,dim_t/])
   ; now load siconc, which we need for masking.
   fin = addfile(b_files_2(n),"r")
   dum_v = fin->$varcode2$
   b_data_2(ind_st:ind_ed,:,:) = dum_v*(conform(dum_v,mask_hemi,(/1,2/)))
   delete([/fin,dum_v/])
end do

; now mask thickness by concentration
b_data = where(b_data_2 .ge. 15.0, b_data, b_data@_FillValue)

;;; bw files
print("Loading bw expt PI control")
bw_data = new((/ntimes_bw,nlat,nlon/),"float")
bw_data_2 = new((/ntimes_bw,nlat,nlon/),"float")
bw_time = new((/ntimes_bw/),"double")

n = 0
do n = 0, nfiles_bw -1
   fin = addfile(bw_files(n),"r")
   ; deal with variable of interest
   dum_v = fin->$varcode$
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
   bw_data(ind_st:ind_ed,:,:) = dum_v
   bw_time(ind_st:ind_ed) = dum_t
   delete([/fin,dum_t,dum_v,dim_t/])
   ; now load siconc, which we need for masking.
   fin = addfile(bw_files_2(n),"r")
   dum_v = fin->$varcode2$
   bw_data_2(ind_st:ind_ed,:,:) = dum_v*(conform(dum_v,mask_hemi,(/1,2/)))
   delete([/fin,dum_v/])
end do

; now mask thickness by concentration
bw_data = where(bw_data_2 .ge. 15.0, bw_data, bw_data@_FillValue)

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
; Get months of interest
;;;;;;;;;;
print ("Getting monthly averages")
b_data_avg = new((/nmons,nlat,nlon/),"float")
bw_data_avg = new((/nmons,nlat,nlon/),"float")

m = 0
do m = 0, nmons -1
   mm = mons(m)
   ; get data
   inds = ind((mm_b .eq. (m+1)).and.((yy_b .ge. 100.0).and.(yy_b .le. 499.0)))
   b_data_avg(m,:,:) = dim_avg_n(b_data(inds,:,:),0)
   delete(inds)
   inds = ind((mm_bw .eq. (m+1)).and.((yy_bw .ge. 100.0).and.(yy_bw .le. 499.0)))
   bw_data_avg(m,:,:) = dim_avg_n(bw_data(inds,:,:),0)
   delete(inds)
   delete(mm)
end do

;;;;;;;;;;
; Get seasons
;;;;;;;;;;
print ("Getting seasonal averages")
nseas = 2
b_data_seas = new((/nseas,nlat,nlon/),"float")
bw_data_seas = new((/nseas,nlat,nlon/),"float")

; b case
;JFM average
b_data_seas(0,:,:) = dim_avg_n(b_data_avg(0:2,:,:),0)
;JAS average
b_data_seas(1,:,:) = dim_avg_n(b_data_avg(6:8,:,:),0)

; bw case
;JFM average
bw_data_seas(0,:,:) = dim_avg_n(bw_data_avg(0:2,:,:),0)
;JAS average
bw_data_seas(1,:,:) = dim_avg_n(bw_data_avg(6:8,:,:),0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Get bins for histogram
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("binning the data")
; make array of levels
nlev = 25
levels = fspan(0,6,nlev)

; now loop through the levels so we can count up how many cells fall
; into each bin. Set to 1 so we can count them up.
; Also need to keep track of areas so we can get fraction of area for
; each thickness level because it's more useful than raw "counts".

b_data_thick_binned = new((/nseas,nlev,nlat,nlon/),typeof(b_data))
b_data_area_binned = new((/nseas,nlev,nlat,nlon/),typeof(b_data))
bw_data_thick_binned = new((/nseas,nlev,nlat,nlon/),typeof(b_data))
bw_data_area_binned = new((/nseas,nlev,nlat,nlon/),typeof(b_data))

; loop through each level
i = 0
do i = 0,nlev-1
   ; how to deal with all levels but last
   if (i .lt. nlev-2)
      ; models
      b_data_thick_binned(:,i,:,:) = where((b_data_seas.gt.levels(i)) .and. (b_data_seas.le.levels(i+1)),1,0)
      b_data_area_binned(:,i,:,:)= b_data_thick_binned(:,i,:,:)*conform(b_data_thick_binned(:,i,:,:),tofloat(tarea_mask_hemi),(/1,2/))
      bw_data_thick_binned(:,i,:,:) = where((bw_data_seas.gt.levels(i)) .and. (bw_data_seas.le.levels(i+1)),1,0)
      bw_data_area_binned(:,i,:,:)= bw_data_thick_binned(:,i,:,:)*conform(bw_data_thick_binned(:,i,:,:),tofloat(tarea_mask_hemi),(/1,2/))
   end if
   ; how to deal with final level
   if (i .eq. nlev-1)
      ; models
      b_data_thick_binned(:,i,:,:) = where((b_data_seas.gt.levels(i)),1,0)
      b_data_area_binned(:,i,:,:)= b_data_thick_binned(:,i,:,:)*conform(b_data_thick_binned(:,i,:,:),tofloat(tarea_mask_hemi),(/1,2/))
      bw_data_thick_binned(:,i,:,:) = where((bw_data_seas.gt.levels(i)),1,0)
      bw_data_area_binned(:,i,:,:)= bw_data_thick_binned(:,i,:,:)*conform(bw_data_thick_binned(:,i,:,:),tofloat(tarea_mask_hemi),(/1,2/))
   end if
end do

print("Calculating area fractions of each binned thicknesses")

; calculate area covered by each ice thickness bin
b_area_binned_summed = dim_sum_n_Wrap(dim_sum_n_Wrap(b_data_area_binned,3),2)
bw_area_binned_summed = dim_sum_n_Wrap(dim_sum_n_Wrap(bw_data_area_binned,3),2)

; calculate total area covered by all ice
b_area_summed = dim_sum_n_Wrap(b_area_binned_summed,1)
bw_area_summed = dim_sum_n_Wrap(bw_area_binned_summed,1)

; calculate fractions
b_spring_area_frac = b_area_binned_summed(0,:)/b_area_summed(0)
b_fall_area_frac = b_area_binned_summed(1,:)/b_area_summed(1)
bw_spring_area_frac = bw_area_binned_summed(0,:)/bw_area_summed(0)
bw_fall_area_frac = bw_area_binned_summed(1,:)/bw_area_summed(1)

; sanity test, these should all sum to 1 (or very very close)
print(sum(b_spring_area_frac))
print(sum(b_fall_area_frac))
print(sum(bw_spring_area_frac))
print(sum(bw_fall_area_frac))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; change levels
levels_plot = levels+0.125  ; need to add (0.25/2) so labels are in the correct spots

fout0 = outdir+"CMIP6_b_bw_"+hemi+"_pi_fractions_for_histogram"

ncdf = addfile(fout0+".nc","c")
ncdf->levels_plot         = levels_plot
ncdf->b_spring_area_frac  = b_spring_area_frac
ncdf->bw_spring_area_frac = bw_spring_area_frac
ncdf->b_fall_area_frac    = b_fall_area_frac
ncdf->bw_fall_area_frac   = bw_fall_area_frac

print("Output file with values")
return

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting PDF of "+varcode)

; set contour levels based on varcode
levels_plot = levels+0.125  ; need to add (0.25/2) so labels are in the correct spots
label_inds = ispan(0,6,1)
label_str = flt2string(int2flt(label_inds))
if (hemi .eq. "NH")
  seas_str = (/"a) Winter - JFM","b) Summer - JAS"/)
end if
if (hemi .eq. "SH")
  seas_str = (/"a) Summer - JFM","b) Winter - JAS"/)
end if

; locations for text
x1 = 5.0
x2 = x1
y1 = 0.3
y2 = 0.27

; set general resources
  res                     = True          ; plot mods desired
  res@gsnDraw             = False         ; don't draw the plot yet
  res@gsnFrame            = False         ; or the frame

; bar chart resources
  res@gsnXYBarChart       = True            ; Create bar plot
  res@gsnXYBarChartOutlineOnly = True 
  res@xyLineThicknessF    = 5.
  res@xyDashPattern       = 0
  res@xyMonoLineColor      = False

; labeling resources
  res@tmXBMode            = "Explicit"         ; explicit labels
  res@tmXBValues          = label_inds
  res@tmXBLabels          = label_str
  res@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  res@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  res@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  res@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis
  res@trXMinF             = 0
  res@trXMaxF             = 6.0
  res@trYMinF             = 0
  res@trYMaxF             = 0.4
  res@tiYAxisString       = "Fraction of total ice area"
  res@tiXAxisString       = "ice thickness bins [m]"

; text resources
  txres = True
  txres@txFontHeightF = 0.03
  txres@txFontThicknessF = 2

; CAM resources
  res_b = res
  res_b@xyLineColors = (/"black"/)
  txres_b = txres
  txres_b@txFontColor = "black"

; WACCM resources
  res_bw = res
  res_bw@xyLineColors = (/"firebrick"/)
  txres_bw = txres
  txres_bw@txFontColor = "firebrick"

;; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@lbOrientation              = "vertical"
  pnlres@gsnPanelFigureStrings      = seas_str
  pnlres@gsnPanelFigureStringsBackgroundFillColor = "white"  ; makes background of label box white
  pnlres@gsnPanelFigureStringsFontHeightF = 0.015   ; text height of label
  pnlres@amJust                     = "TopRight"                 ; put label in top right corner
  pnlres@gsnPanelFigureStringsPerimOn = False ; no box around labels

;;;;;;;;;;;;;
; Make regular plots
;;;;;;;;;;;;;
print("Plotting seasonal histograms")
; set info for this plot
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_histogram"

  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,"wh-bl-gr-ye-re") 

; do individual plots
plot = new(2,graphic)
s = 0
bplot = gsn_csm_xy(wks,levels_plot,b_spring_area_frac,res_b)
oplot1 = gsn_csm_xy(wks,levels_plot,bw_spring_area_frac,res_bw)
overlay(bplot,oplot1)
dum1 = gsn_add_text(wks,bplot,"CAM",x1,y1,txres_b)
dum2 = gsn_add_text(wks,bplot,"WACCM",x2,y2,txres_bw)
plot(s) = bplot

s = 1
bplot = gsn_csm_xy(wks,levels_plot,b_fall_area_frac,res_b)
oplot1 = gsn_csm_xy(wks,levels_plot,bw_fall_area_frac,res_bw)
overlay(bplot,oplot1)
plot(s) = bplot

; add text overlay

; panel plot
print("Paneling")
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,2/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Done with all plots. Hooray!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
