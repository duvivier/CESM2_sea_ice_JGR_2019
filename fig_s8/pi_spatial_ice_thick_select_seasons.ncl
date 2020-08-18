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
vartitle = "ice_thick"
ti_str = "PI-Control: Ice Thickness"
ti_str2 = ti_str+"~C~Std.dev"
ti_str3 = ti_str+"~C~Diff(CAM-WACCM)"
units = "m"

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

mons = fspan(1.,12.,12)
nmons = dimsizes(mons)
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
   ;fin = addfile(data_dir+"REGION_MASK_gx1v6.nc","r")
   ;mask_in = fin->ArcOc_mask(0,:,:)
   fin = addfile(data_dir+"CCSM4_Arctic_mask.nc","r")
   mask_in = fin->CCSM4_Arctic_mask_no_Barents_Kara_CAA_Nares
   delete(fin)
   ; whole hemisphere
   ;mask_in = where(lat2d .ge. 70.0, 1.0, 0.0)
end if
if (hemi .eq. "SH")then
   mask_in = where(lat2d .le. -50.0, 1.0, 0.0)
end if

; reset missing value info
mask_in@_FillValue = -9999
mask_in = where(ismissing(mask_in), 0.0, mask_in)

; get hemispheric mask of areas
; units are m^2
mask_hemi = tofloat(tarea)*mask_in

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

; now mask thickness by concentration
b_data = where(b_data_2_all .ge. 15.0, b_data_all, b_data_all@_FillValue)

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

; now mask thickness by concentration
bw_data = where(bw_data_2_all .ge. 15.0, bw_data_all, bw_data_all@_FillValue)

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
; Get Seasonal averages
;;;;;;;;;;
print ("Getting seasonal averages")
nseas = 4
b_data_seas = new((/nseas,nlat,nlon/),"float")
bw_data_seas = new((/nseas,nlat,nlon/),"float")
b_data_seas_var = new((/nseas,nlat,nlon/),"float")
bw_data_seas_var = new((/nseas,nlat,nlon/),"float")

s = 0
do s = 0, nseas -1
  ; b-case data
  ; JFM
  if (s .eq. 0)
     inds = ind(((mm_b .eq. 1) .or. (mm_b .eq. 2) .or. (mm_b .eq. 3)) .and. ((yy_b .ge. 100.0).and.(yy_b .le.499)))
  end if
  ; AMJ
  if (s .eq. 1)
     inds = ind(((mm_b .eq. 4) .or. (mm_b .eq. 5) .or. (mm_b .eq. 6)) .and. ((yy_b .ge. 100.0).and.(yy_b .le.499)))
  end if
  ; JAS
  if (s .eq. 2)
     inds = ind(((mm_b .eq. 7) .or. (mm_b .eq. 8) .or. (mm_b .eq. 9)) .and. ((yy_b .ge. 100.0).and.(yy_b .le.499)))
  end if
  ; OND
  if (s .eq. 3)
     inds = ind(((mm_b .eq. 10) .or. (mm_b .eq. 11) .or. (mm_b .eq. 12)) .and. ((yy_b .ge. 100.0).and.(yy_b .le.499)))
  end if
     ntot_b = dimsizes(inds)
     dum = b_data(inds,:,:) ;reshape(b_data(inds,:,:),(/ntot_b,nlat,nlon/))
     b_data_seas(s,:,:) = dim_avg_n_Wrap(dum,0)
     b_data_seas_var(s,:,:) = dim_variance_n_Wrap(dum,0)
     delete([/inds,dum/])
  ; bw-case data
  ; JFM
  if (s .eq. 0)
     inds = ind(((mm_bw .eq. 1) .or. (mm_bw .eq. 2) .or. (mm_bw .eq. 3)) .and. ((yy_bw .ge. 100.0).and.(yy_bw .le.499)))
  end if
  ; AMJ
  if (s .eq. 1)
     inds = ind(((mm_bw .eq. 4) .or. (mm_bw .eq. 5) .or. (mm_bw .eq. 6)) .and. ((yy_bw .ge. 100.0).and.(yy_bw .le.499)))
  end if
  ; JAS
  if (s .eq. 2)
     inds = ind(((mm_bw .eq. 7) .or. (mm_bw .eq. 8) .or. (mm_bw .eq. 9)) .and. ((yy_bw .ge. 100.0).and.(yy_bw .le.499)))
  end if
  ; OND
  if (s .eq. 3)
     inds = ind(((mm_bw .eq. 10) .or. (mm_bw .eq. 11) .or. (mm_bw .eq. 12)) .and. ((yy_bw .ge. 100.0).and.(yy_bw .le.499)))
  end if
     ntot_bw = dimsizes(inds)
     dum = bw_data(inds,:,:) ;reshape(b_data(inds,:,:),(/ntot_bw,nlat,nlon/))
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
; Get differences
;;;;;;;;;;
if (varcode .eq. "sivol")
   n = 0
   do n = 0, nseas -1
      b_data_seas(n,:,:) = where(((kmt .gt. 0.5).and.(ismissing(b_data_seas(n,:,:)))), -0.0000001, b_data_seas(n,:,:))
      bw_data_seas(n,:,:) = where(((kmt .gt. 0.5).and.(ismissing(bw_data_seas(n,:,:)))), -0.0000001, bw_data_seas(n,:,:))
   end do
end if

; calculate difference
diff_seas = b_data_seas - bw_data_seas

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set contour levels based on varcode
if (varcode .eq. "sivol")
   main_cntrs=(/0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0/)
   diff_cntrs=(/-1.5,-1.0,-0.75,-0.5,-0.25,-0.1,0.1,0.25,0.5,0.75,1.0,1.5/)
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
  res@cnMissingValFillColor = "grey"
  res@cnLevelSelectionMode    = "ExplicitLevels"

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
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_seasonal"
   title_fin = ti_str+ " ("+units+") WACCM"

; open wks
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")

; set resources
  res_main                = res
  res_main@cnLevels       = main_cntrs

; do individual plots
plot = new(2, graphic)
; Winter
s = 0
res_main@gsnCenterString = seas_str(s)
plot(s) = gsn_csm_contour_map(wks,bw_data_seas(0,:,:),res_main)
; Summer
s = 1
res_main@gsnCenterString = seas_str(s)
plot(s) = gsn_csm_contour_map(wks,bw_data_seas(2,:,:),res_main)

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
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_diff_seasonal"
   title_fin = ti_str + " ("+units+") difference: CAM-WACCM"

; open wks
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
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
