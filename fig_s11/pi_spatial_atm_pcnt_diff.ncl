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
;varcode = "clt"
; "clt" "clivi"
;hemi = "NH"
; ******************************************************
; BEGIN SCRIPT
;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
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
nlat = dimsizes(lat1d)
nlon = dimsizes(lon1d)
delete([/b_file,dummy2d/])

mons = fspan(1.,12.,12)
nmons = dimsizes(mons)
nyrs_b = 1200
nyrs_bw = 499
ntimes_b = nmons*nyrs_b
ntimes_bw = nmons*nyrs_bw

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data")
; data is in 50 or 100 year chunks

;;; b files
print("Loading b expt PI control")
b_data = new((/ntimes_b,nlat,nlon/),"float")
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
end do

;;; bw files
print("Loading bw expt PI control")
bw_data = new((/ntimes_bw,nlat,nlon/),"float")
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
end do

; convert Pa to hPa for slp
if (varcode .eq. "psl")
   b_data = b_data/100.
   bw_data = bw_data/100.
   b_data@units = "hPa"
end if

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
   inds = ind((mm_b .eq. (m+1)).and.((yy_b .ge. 100.0).and.(yy_b.le. 499.0)))
   b_data_avg(m,:,:) = dim_avg_n(b_data(inds,:,:),0)
   delete(inds)
   inds = ind((mm_bw .eq. (m+1)).and.((yy_bw .ge. 100.0).and.(yy_bw.le. 499.0)))
   bw_data_avg(m,:,:) = dim_avg_n(bw_data(inds,:,:),0)
   delete(inds)   
   delete(mm)
end do

; calculate percent difference
diff_avg = ((b_data_avg - bw_data_avg)/bw_data_avg)*100.

;;;;;;;;;;
; Assign lat/lon for plotting
;;;;;;;;;;
diff_avg!0 = "mons"
diff_avg&mons = mons
diff_avg!1 = "lat"
diff_avg&lat = lat1d
diff_avg!2 = "lon"
diff_avg&lon = lon1d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting "+varcode)

; set contour levels based on varcode
if (varcode .eq. "clivi")
   diff_cntrs = (/-40.,-30.,-20.,-10.,0.,10.,20.,30.,40./)
end if
if (varcode .eq. "clt")
   diff_cntrs = (/-8.,-6.,-4.,-2.,0.,2.,4.,6.,8./)
end if
mons_str = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct","Nov","Dec"/)

; set outdir and fout names
   outdir = "./"
   fout = outdir+"CMIP6_b_bw_"+hemi+"_piControl_"+varcode+"_pcnt_diff_monthly"
   title_fin = b_data@title + " ("+b_data@units+") difference: (CAM-WACCM)/WACCM"

; output file
  wks_type = "ps" 
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,"GMT_polar")

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
  res@mpGridAndLimbOn      = False    ; set to true to show lat and lon lines
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
  res@cnMissingValFillColor = "white"
  res@cnLevelSelectionMode    = "ExplicitLevels"

; set diff contour resources
  resdiff                     = res
  resdiff@cnLevels            = diff_cntrs

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

;;;;;;;;;;;;;
; Make individual plots
;;;;;;;;;;;;;
print ("Plotting "+varcode)
plot = new(5, graphic)
m = 0
do m = 0, 4
   if (m .eq. 0)
      mm = 4
   end if
   resdiff@gsnCenterString = mons_str(mm)
   plot(m) = gsn_csm_contour_map(wks,diff_avg(mm,:,:),resdiff)
   mm = mm + 1
end do

;;;;;;;;;;;;;
; Panel plot
;;;;;;;;;;;;;
print("Paneling")
  ; make panel plot
  ;pnlres@txString  = title_fin
  gsn_panel(wks,(/plot/),(/1,5/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Done with all plots. Hooray!")

;;;;;;;;;;;;;;;;;;;;;; END script
end
