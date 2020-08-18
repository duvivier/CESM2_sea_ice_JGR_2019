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
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
;***************************************
begin
;;;;;;;;;;
; Manual inputs or from cshell script
;;;;;;;;;;
; ******************************************************
; BEGIN SCRIPT

; load in all the variables
varcode1  = "FSNT"     ; net SW @TOM 
varcode2  = "FLNT"     ; net LW @TOM
varcode3  = "FSNS"     ; net SW @SFC
varcode4  = "FLNS"     ; net LW @SFC
varcode5  = "SHFLX"    ; sensible heat flux @SFC
varcode6  = "LHFLX"    ; latent heat flux @SFC (needs correction)
varcode7  = "PRECSC"   ; convective snow precip
varcode8  = "PRECSL"   ; large scale snow precip
varcode9  = "N_HEAT"   ; ocean diagnostic heat transport
varcode10 = "dvidtd"   ; volume tendency due to dynamics

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
cmip_dir = "/glade/collections/cdg/timeseries-cmip6/"
grid_dir = "/glade/p/cesm/omwg/grids/"

; case 1
case1 = "b.e21.B1850.f09_g17.CMIP6-piControl.001"
b_files1 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*."+varcode1+".*.nc")
b_files2 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*."+varcode2+".*.nc")
b_files3 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*."+varcode3+".*.nc")
b_files4 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*."+varcode4+".*.nc")
b_files5 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*."+varcode5+".*.nc")
b_files6 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*."+varcode6+".*.nc")
b_files7 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*."+varcode7+".*.nc")
b_files8 = systemfunc("ls "+cmip_dir+case1+"/atm/proc/tseries/month_1/*."+varcode8+".*.nc")
b_files9 = systemfunc("ls "+cmip_dir+case1+"/ocn/proc/tseries/month_1/*."+varcode9+".*.nc")
b_files10 = systemfunc("ls "+cmip_dir+case1+"/ice/proc/tseries/month_1/*."+varcode10+".*.nc")
nfiles_b = dimsizes(b_files1)

; case 2
case2 = "b.e21.BW1850.f09_g17.CMIP6-piControl.001"
bw_files1 = systemfunc("ls "+cmip_dir+case2+"/atm/proc/tseries/month_1/*."+varcode1+".*.nc")
bw_files2 = systemfunc("ls "+cmip_dir+case2+"/atm/proc/tseries/month_1/*."+varcode2+".*.nc")
bw_files3 = systemfunc("ls "+cmip_dir+case2+"/atm/proc/tseries/month_1/*."+varcode3+".*.nc")
bw_files4 = systemfunc("ls "+cmip_dir+case2+"/atm/proc/tseries/month_1/*."+varcode4+".*.nc")
bw_files5 = systemfunc("ls "+cmip_dir+case2+"/atm/proc/tseries/month_1/*."+varcode5+".*.nc")
bw_files6 = systemfunc("ls "+cmip_dir+case2+"/atm/proc/tseries/month_1/*."+varcode6+".*.nc")
bw_files7 = systemfunc("ls "+cmip_dir+case2+"/atm/proc/tseries/month_1/*."+varcode7+".*.nc")
bw_files8 = systemfunc("ls "+cmip_dir+case2+"/atm/proc/tseries/month_1/*."+varcode8+".*.nc")
bw_files9 = systemfunc("ls "+cmip_dir+case2+"/ocn/proc/tseries/month_1/*."+varcode9+".*.nc")
bw_files10 = systemfunc("ls "+cmip_dir+case2+"/ice/proc/tseries/month_1/*."+varcode10+".*.nc")
nfiles_bw = dimsizes(bw_files1)

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get atmosphere dimensions
fin = addfile(b_files1(0),"r")
lat1d_atm = fin->lat
lon1d_atm = fin->lon
gw = tofloat(fin->gw)
nlat_atm = dimsizes(lat1d_atm)
nlon_atm = dimsizes(lon1d_atm)
delete(fin)

; get ocean dimensions
fin = addfile(grid_dir+"gx1v7_grid.nc","r")
lat2d_ocn = fin->TLAT
lon2d_ocn = fin->TLONG
dims = dimsizes(lat2d_ocn)
nlat_ocn = dims(0)
nlon_ocn = dims(1)
delete(fin)
; get some dimensions for N_HEAT file
nreg = 2
ncomp = 5
nlat_ocn_aux = 395

nmons = 12
nyrs_b = 1200
nyrs_bw = 499
ntimes_b = nmons*nyrs_b
ntimes_bw = nmons*nyrs_bw

;;;;;;;;;;
; Set constants
;;;;;;;;;;
    pi      = 3.14159265    ; pi
    g       = 9.81          ; gravity (m/s2)
    Lf      = 3.337e5 	    ; latent heat of fusion (J/kg)
    Re      = 6.371e6       ; radius of earth (m)
    rho_h20 = 1000.         ; density of water (kg/m3)
    rho_ice = 917.          ; density of ice (kg/m3)

    ; derived values for NHT
    dlon = (2.*pi)/nlon_atm          ; change in longitude (radians) per unit latitude
    coef = (Re^2)/1.e15             ; (m2 and scaled for PW)
    ; derived values for ice NHT
    Earea = 4.*pi*Re^2               ; surface area of Earth (m2)
    dlon_deg = (lon1d_atm(2) - lon1d_atm(1))  ; assume dlon is constant
    dlat_deg = (lat1d_atm(2) - lon1d_atm(1))  ; assume dlat is constant
    wgt = tofloat(NormCosWgtGlobe(lat1d_atm))
    tsum = sum(wgt)          ; sum of all weights
    nwgt = wgt/tsum          ; fraction of sphre of each latitude band
    boxarea = Earea*nwgt/nlon_atm  ; area of each grid box latitude in m2

    ; get surface area of Earth above each latitude band
    ; A = 2*pi*Re^2*(1-cos(90-lat))  - NOTE, must convert latitude to radians
    lat_in = tofloat((90.0-lat1d_atm)*(pi/180.))  ; get latitude to input in radians
    area_sfc_cap = (2*pi*(Re^2))*(1-cos(lat_in))  ; in m2

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
print("Loading CESM data, Calculate hemispheric means at the same time")
; data is in 50 or 100 year chunks
; avg over just hemi cells

;;; b files
print("Loading b expt PI control")
b_data_1 = new((/ntimes_b,nlat_atm,nlon_atm/),"float")
b_data_2 = new((/ntimes_b,nlat_atm,nlon_atm/),"float")
b_data_3 = new((/ntimes_b,nlat_atm,nlon_atm/),"float")
b_data_4 = new((/ntimes_b,nlat_atm,nlon_atm/),"float")
b_data_5 = new((/ntimes_b,nlat_atm,nlon_atm/),"float")
b_data_6 = new((/ntimes_b,nlat_atm,nlon_atm/),"float")
b_data_7 = new((/ntimes_b,nlat_atm,nlon_atm/),"float")
b_data_8 = new((/ntimes_b,nlat_atm,nlon_atm/),"float")
b_data_9 = new((/ntimes_b,nreg,ncomp,nlat_ocn_aux/),"float")
b_data_10 = new((/ntimes_b,nlat_ocn,nlon_ocn/),"float")
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
   b_data_1(ind_st:ind_ed,:,:) = fin1->$varcode1$
   fin2 = addfile(b_files2(n),"r")
   b_data_2(ind_st:ind_ed,:,:) = fin2->$varcode2$
   fin3 = addfile(b_files3(n),"r")
   b_data_3(ind_st:ind_ed,:,:) = fin3->$varcode3$
   fin4 = addfile(b_files4(n),"r")
   b_data_4(ind_st:ind_ed,:,:) = fin4->$varcode4$
   fin5 = addfile(b_files5(n),"r")
   b_data_5(ind_st:ind_ed,:,:) = fin5->$varcode5$
   fin6 = addfile(b_files6(n),"r")
   b_data_6(ind_st:ind_ed,:,:) = fin6->$varcode6$
   fin7 = addfile(b_files7(n),"r")
   b_data_7(ind_st:ind_ed,:,:) = fin7->$varcode7$
   fin8 = addfile(b_files8(n),"r")
   b_data_8(ind_st:ind_ed,:,:) = fin8->$varcode8$
   fin9 = addfile(b_files9(n),"r")
   b_data_9(ind_st:ind_ed,:,:,:) = fin9->$varcode9$
   fin10 = addfile(b_files10(n),"r")
   b_data_10(ind_st:ind_ed,:,:) = fin10->$varcode10$
   delete([/fin1,fin2,fin3,fin4,fin5,fin6,fin7,fin8,fin9,fin10,dum_t/])
end do

;;; bw files
print("Loading bw expt PI control")
bw_data_1 = new((/ntimes_bw,nlat_atm,nlon_atm/),"float")
bw_data_2 = new((/ntimes_bw,nlat_atm,nlon_atm/),"float")
bw_data_3 = new((/ntimes_bw,nlat_atm,nlon_atm/),"float")
bw_data_4 = new((/ntimes_bw,nlat_atm,nlon_atm/),"float")
bw_data_5 = new((/ntimes_bw,nlat_atm,nlon_atm/),"float")
bw_data_6 = new((/ntimes_bw,nlat_atm,nlon_atm/),"float")
bw_data_7 = new((/ntimes_bw,nlat_atm,nlon_atm/),"float")
bw_data_8 = new((/ntimes_bw,nlat_atm,nlon_atm/),"float")
bw_data_9 = new((/ntimes_bw,nreg,ncomp,nlat_ocn_aux/),"float")
bw_data_10 = new((/ntimes_bw,nlat_ocn,nlon_ocn/),"float")
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
   bw_data_1(ind_st:ind_ed,:,:) = fin1->$varcode1$
   fin2 = addfile(bw_files2(n),"r")
   bw_data_2(ind_st:ind_ed,:,:) = fin2->$varcode2$
   fin3 = addfile(bw_files3(n),"r")
   bw_data_3(ind_st:ind_ed,:,:) = fin3->$varcode3$
   fin4 = addfile(bw_files4(n),"r")
   bw_data_4(ind_st:ind_ed,:,:) = fin4->$varcode4$
   fin5 = addfile(bw_files5(n),"r")
   bw_data_5(ind_st:ind_ed,:,:) = fin5->$varcode5$
   fin6 = addfile(bw_files6(n),"r")
   bw_data_6(ind_st:ind_ed,:,:) = fin6->$varcode6$
   fin7 = addfile(bw_files7(n),"r")
   bw_data_7(ind_st:ind_ed,:,:) = fin7->$varcode7$
   fin8 = addfile(bw_files8(n),"r")
   bw_data_8(ind_st:ind_ed,:,:) = fin8->$varcode8$
   fin9 = addfile(bw_files9(n),"r")
   bw_data_9(ind_st:ind_ed,:,:,:) = fin9->$varcode9$
   fin10 = addfile(bw_files10(n),"r")
   bw_data_10(ind_st:ind_ed,:,:) = fin10->$varcode10$
   delete([/fin1,fin2,fin3,fin4,fin5,fin6,fin7,fin8,fin9,fin10,dum_t/])
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

; get just indices of overlapping time
inds_b = ind((yy_b .ge. 100.0).and.(yy_b .le. 499.0))
inds_bw = ind((yy_bw .ge. 100.0).and.(yy_bw .le. 499.0))

;;;;;;;;;;
; Calculate total NHT
;;;;;;;;;;
print("Finding total Northward Heat Transport")
; Uses implied total NW heat transport
; Use method from Kay et al. 2012, J.Clim

; Get top of model net radiative fluxes
b_net_toa = b_data_1(inds_b,:,:) - b_data_2(inds_b,:,:)
bw_net_toa = bw_data_1(inds_bw,:,:) - bw_data_2(inds_bw,:,:)

; average over times, sum over longitudes
b_heatflux = dim_sum_n(dim_avg_n(b_net_toa,0),1)
bw_heatflux = dim_sum_n(dim_avg_n(bw_net_toa,0),1)

; now loop through latitudes to find value at each latitude
; this is the total NHT required by top of model values
; units are PetaWatts ; dimension is lat
b_nht_tot = new((/nlat_atm/),"float")
bw_nht_tot = new((/nlat_atm/),"float")

; start at northernmost point
do j = nlat_atm-1, 0, 1
   b_nht_tot(j) = -coef*dlon*sum(b_heatflux(j:nlat_atm-1)*gw(j:nlat_atm-1))
   bw_nht_tot(j) = -coef*dlon*sum(bw_heatflux(j:nlat_atm-1)*gw(j:nlat_atm-1))
end do

delete([/b_heatflux,bw_heatflux/])

;;;;;;;;;;
; Calculate atmosphere NHT
;;;;;;;;;;
print("Finding atmospheric Northward Heat Transport")
; Uses implied total NW heat transport
; Use method from Kay et al. 2012, J.Clim

; first, need to correct latent heat flux to account for snow
; use method developed by kay/holland. Units are W/m2
b_data_6_corr = b_data_6(inds_b,:,:) + (Lf*rho_h20*(b_data_7(inds_b,:,:)+b_data_8(inds_b,:,:)))
bw_data_6_corr = bw_data_6(inds_bw,:,:) + (Lf*rho_h20*(bw_data_7(inds_bw,:,:)+bw_data_8(inds_bw,:,:)))

; Get surface net fluxes
; Note that positive for LW, SH, and LH mean into ATM
b_net_sfc = b_data_3(inds_b,:,:) - b_data_4(inds_b,:,:) - b_data_5(inds_b,:,:) - b_data_6_corr
bw_net_sfc = bw_data_3(inds_bw,:,:) - bw_data_4(inds_bw,:,:) - bw_data_5(inds_bw,:,:) - bw_data_6_corr

; get difference between TOM and SFC
b_net = (-b_net_toa + b_net_sfc)
bw_net = (-bw_net_toa + bw_net_sfc)

; average over times, sum over longitudes
b_heatflux = dim_sum_n(dim_avg_n(b_net,0),1)
bw_heatflux = dim_sum_n(dim_avg_n(bw_net,0),1)

; now loop through latitudes to find value at each latitude
; this is the atmospheric component of the NHT
; units are PetaWatts ; dimension is lat
b_nht_atm = new((/nlat_atm/),"float")
bw_nht_atm = new((/nlat_atm/),"float")

; start at northernmost point
do j = nlat_atm-1, 0, 1
   b_nht_atm(j) = coef*dlon*sum(b_heatflux(j:nlat_atm-1)*gw(j:nlat_atm-1))
   bw_nht_atm(j) = coef*dlon*sum(bw_heatflux(j:nlat_atm-1)*gw(j:nlat_atm-1))
end do

delete([/b_heatflux,bw_heatflux/])

;;;;;;;;;;
; Calculate ocean NHT
;;;;;;;;;;
print("Finding oceanic Northward Heat Transport")
; Uses pop diagnostic NHT output

; get dimensions
; data should be: time x transport_region x transport_component x latitudes
lat1d_ocn = b_data_9&lat_aux_grid
;regions = b_data_9&transport_regions  ;want region 0 - Global Oceans and Marginal Seas
;components = b_data_9&transport_components  ;want component 0 - total

; average over time for relevant subset (global, total)
b_data_9_sub = dim_avg_n(b_data_9(inds_b,0,0,:),0)
bw_data_9_sub = dim_avg_n(bw_data_9(inds_bw,0,0,:),0)

; now convert to the atmosphere grid
b_nht_ocn = linint1(lat1d_ocn,b_data_9_sub,False,lat1d_atm,0)
bw_nht_ocn = linint1(lat1d_ocn,bw_data_9_sub,False,lat1d_atm,0)

;;;;;;;;;;
; Calculate ice NHT
;;;;;;;;;;
print("Finding sea ice Northward Heat Transport")
;; use the volume convergence (volume tendency) due to dynamics
;; to calculate heat transport from ice export out of Arctic ocean

;; average over time
;b_data_10_sub = dim_avg_n(b_data_10(inds_b,:,:),0)
;bw_data_10_sub = dim_avg_n(bw_data_10(inds_bw,:,:),0)

;; Interpolate to the atmosphere grid
;b_data_10_atm = PopLatLon(b_data_10_sub,"gx1v6","fv0.9x1.25","aave","da","090309")
;bw_data_10_atm = PopLatLon(bw_data_10_sub,"gx1v6","fv0.9x1.25","aave","da","090309")

;; Sum over the longitudes
;b_data_10_atm_sub = dim_sum_n_Wrap(b_data_10_atm,1)
;bw_data_10_atm_sub = dim_sum_n_Wrap(bw_data_10_atm,1)

;; multiply by area of grid cell to convert to ice flux of cm/day*m2
;b_data_10_atm_sub = b_data_10_atm_sub/boxarea
;bw_data_10_atm_sub = bw_data_10_atm_sub/boxarea
;; convert to m3/s --> 1m/100cm,1day/24hr,1hr/60m,60s/1m
;b_data_10_atm_sub = b_data_10_atm_sub/(100.*24.*60.*60.)
;bw_data_10_atm_sub = bw_data_10_atm_sub/(100.*24.*60.*60.)

;; find northward ice flux
;b_north_ice_flux = b_data_10_atm_sub*0.0
;bw_north_ice_flux = bw_data_10_atm_sub*0.0
;j = 0
;do j = 0,nlat_atm-2
;   b_north_ice_flux(j) = b_data_10_atm_sub(j) + sum(b_data_10_atm_sub(j+1:nlat_atm-1))
;   bw_north_ice_flux(j) = bw_data_10_atm_sub(j) + sum(bw_data_10_atm_sub(j+1:nlat_atm-1))
;end do

;; convert from northward ice flux (m^3/sec) to heat flux(J/s = W)
;b_nht_ice = -(b_north_ice_flux*rho_ice*Lf*to_PW)
;bw_nht_ice = -(bw_north_ice_flux*rho_ice*Lf*to_PW)

; calculate as a residual for now
b_nht_ice = b_nht_tot-(b_nht_atm+b_nht_ocn)
bw_nht_ice = bw_nht_tot-(bw_nht_atm+bw_nht_ocn)

;;;;;;;;;;
; Combine for plotting
;;;;;;;;;;
b_data_all = new((/5,nlat_atm/),typeof(b_data_1))
bw_data_all = new((/5,nlat_atm/),typeof(bw_data_1))

b_data_all(0,:) = b_nht_tot     ; total
b_data_all(1,:) = b_nht_atm     ; atm
b_data_all(2,:) = b_nht_ocn     ; ocn
b_data_all(3,:) = b_nht_ice     ; implied ice
;b_data_all(4,:) = b_nht_atm+b_nht_ocn

bw_data_all(0,:) = bw_nht_tot
bw_data_all(1,:) = bw_nht_atm
bw_data_all(2,:) = bw_nht_ocn
bw_data_all(3,:) = bw_nht_ice
;bw_data_all(4,:) = bw_nht_atm+bw_nht_ocn

;;;;;;;;;;
; Get differences
;;;;;;;;;;
diff = b_data_all - bw_data_all
pcnt_diff = ((b_data_all - bw_data_all)/bw_data_all)*100.
pcnt_diff = where((diff .lt. 0.001),0.0, pcnt_diff)

;;;;;;;;;;
; Also get in W/m2
;;;;;;;;;;
; convert to W/m2 instead of PW so we can directly compare to the radiative forcing
; We want to divide the NHT by the area north of a given latitude since the
; heat is distributed "evenly" over that whole area
b_data_all_wm2 = new(dimsizes(b_data_all),typeof(b_data_all))
bw_data_all_wm2 = new(dimsizes(bw_data_all),typeof(bw_data_all))

; how to deal with near-zero values
area_sfc_cap = where(area_sfc_cap .eq. 0.0, area_sfc_cap@_FillValue, area_sfc_cap)
n = 0
do n = 0, 4
   b_data_all_wm2(n,:) = (b_data_all(n,:)*1.0e15)/area_sfc_cap
   bw_data_all_wm2(n,:) = (bw_data_all(n,:)*1.0e15)/area_sfc_cap
end do

diff_wm2 = b_data_all_wm2 - bw_data_all_wm2
pcnt_diff_wm2 = ((b_data_all_wm2 - bw_data_all_wm2)/bw_data_all_wm2)*100.
pcnt_diff_wm2 = where((diff_wm2 .lt. 0.001),0.0, pcnt_diff_wm2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together with other data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
outdir = "./"
fout0 = outdir+"CMIP6_b_bw_northward_heat_transport"

ncdf = addfile(fout0+".nc","c")
ncdf->b_data_all          = b_data_all
ncdf->bw_data_all         = bw_data_all
ncdf->diff                = diff
ncdf->b_data_all_wm2      = b_data_all_wm2
ncdf->bw_data_all_wm2     = bw_data_all_wm2
ncdf->diff_wm2            = diff_wm2 

return
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting NHT")
; set outdir and fout names
   outdir = "./"

units = "PW"

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
  opts@trXMinF                     = -90.0
  opts@trXMaxF                     = 90.0
  opts@tiXAxisString               = "Latitude"
  opts@tiYAxisString               = units
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

; resources for b
  opts_b = opts
  opts_b@tiMainOn = False
  opts_b@xyMarkLineMode = "Lines"
  opts_b@xyLineColors = (/"black","firebrick","blue","deepskyblue"/)
  opts_b@xyLineThicknesses = (/5.0,4.0,4.0,4.0/)
  opts_b@xyDashPatterns = (/0,0,0,0/)

; text
  title1 = "Net NHT"
  txres_1 = txres
  txres_1@txFontColor = "black"
  title2 = "ATM NHT"
  txres_2 = txres
  txres_2@txFontColor = "firebrick"
  title3 = "OCN NHT"
  txres_3 = txres
  txres_3@txFontColor = "blue"
  title4 = "ICE NHT"
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
; Make b case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_full_piControl_NHT"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)  

   ymin = -6.0
   ymax = 6.0
   yminor = 3
   x0 = -65.0
   x1 = -60.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 5.5
   y1 = 4.0
   y2 = 3.25
   y3 = 2.5
   y4 = 1.75 

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor  

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,b_data_all,opts_b)
   dum0 = gsn_add_text(wks,lplot,"a) NHT CAM6",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
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

;;;;;;;;;;;;;
; Make bw case plot
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_bw_full_piControl_NHT"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)    

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,bw_data_all,opts_b)
   dum0 = gsn_add_text(wks,lplot,"a) NHT WACCM6",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
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

;;;;;;;;;;;;;
; Make difference plots
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_full_piControl_NHT_diffs"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)

   ymin = -0.1
   ymax = 0.1
   yminor = 3
   x0 = 0.0 ;-60.0
   x1 = -40.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 0.085
   y1 = -0.04
   y2 = -0.05
   y3 = -0.06
   y4 = -0.07

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor   

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,diff,opts_b)
   ;dum0 = gsn_add_text(wks,lplot,"b) NHT Difference~C~(CAM6-WACCM6)",x0,y0,txres_1)
   dum0 = gsn_add_text(wks,lplot,"b) NHT Difference",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
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

;;;;;;;;;;;;;
; Make percent difference plots
;;;;;;;;;;;;;
 fout = outdir+"CMIP6_b_bw_full_piControl_NHT_pcnt_diffs"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)  

   ymin = -15.0
   ymax = 15.0
   yminor = 4
   x0 = 65.0
   x1 = 60.0
   x2 = x1
   x3 = x1
   x4 = x1
   y0 = 14.0
   y1 = 10.0
   y2 = 8.0
   y3 = 6.0
   y4 = 4.0

  opts_b@trYMinF                     = ymin
  opts_b@trYMaxF                     = ymax
  opts_b@tmYLMinorPerMajor           = yminor
  opts_b@tmYRMinorPerMajor           = yminor    
  opts@tiYAxisString                 = "% diff"

   ; do actual plotting
   lplot = gsn_csm_xy(wks,lat1d_atm,pcnt_diff,opts_b)
   dum0 = gsn_add_text(wks,lplot,"(CAM-WACCM)/WACCM",x0,y0,txres_1)
   dum1 = gsn_add_text(wks,lplot,title1,x1,y1,txres_1)
   dum2 = gsn_add_text(wks,lplot,title2,x2,y2,txres_2)
   dum3 = gsn_add_text(wks,lplot,title3,x3,y3,txres_3)
   dum4 = gsn_add_text(wks,lplot,title4,x4,y4,txres_4)
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
