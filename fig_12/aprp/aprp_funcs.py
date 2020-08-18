#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Functions for partitioning shortwave using APRP (Taylor et al. 2007)
#Coding by Jennifer Kay (Jennifer.E.Kay@colorado.edu) and Elizabeth Maroon (emaroon@wisc.edu)
#Last updated: November 13, 2019

import numpy as np
import xarray as xr
import glob
import matplotlib.pyplot as plt

def read_data(cnt_name,cnt_path,cnt_startyr,cnt_startyr_APRP, \
                      exp_name,exp_path,exp_startyr,exp_startyr_APRP,Nyrs_APRP,Ftype):

    ## define variables to load and files that contain them -- Load them :)
    ds_list={}
    #vars_to_load=['SOLIN'] ## use when testing inputs
    vars_to_load=['SOLIN','FSDS','FSDSC','FSNT','FSNS','FSNTC','FSNSC','CLDTOT','TS']
    expnames=[]
    expnames.append('control')
    expnames.append('experiment')
    expfiles={}

    ds_cnt={}
    ds_exp={}

    for vv in vars_to_load:
        print("working on "+vv)

        #### Find and check file names
        if Ftype=='le':
            expfiles[expnames[0]]=cnt_path+vv+'/'+cnt_name+'.cam.h0.'
            #print(expfiles[expnames[0]])
            filename_cnt=glob.glob(expfiles[expnames[0]]+vv+'.'+cnt_startyr+'*.nc')[0]
            #print(filename_cnt)
        if Ftype=='default':    
            expfiles[expnames[0]]=cnt_path+'/'+cnt_name+'.cam.h0.'
            #print(expfiles[expnames[0]])
            filename_cnt=glob.glob(expfiles[expnames[0]]+vv+'.'+cnt_startyr+'*.nc')[0]
            #print(filename_cnt)        
        expfiles[expnames[1]]=exp_path+'/'+exp_name+'.cam.h0.'
        #print(expfiles[expnames[1]])
        #load in experiment arrays
        filename_exp=glob.glob(expfiles[expnames[1]]+vv+'.'+exp_startyr+'*.nc')[0]
        #print(filename_exp)

        ### Load in control (ds_cnt)
        if Ftype=='le':
            #print(filename_cnt)
            ds_temp=xr.open_dataset(filename_cnt,decode_times=False)
            #print(ds_temp['time'][0],ds_temp['time'][-1])
            ## what is this 15? I think it is a 15 day rewind for CAM history files.  Double check with Elizabeth.
            ds_temp['time'].values=ds_temp['time'].values-15 
            ds_temp=xr.decode_cf(ds_temp,decode_times=True)
            #print(ds_temp['time'][0],ds_temp['time'][-1])    
            numyrs=int(np.floor(len(ds_temp['time'])/12))
            cnt_endyr_APRP=(cnt_startyr_APRP+Nyrs_APRP-int(cnt_startyr))
            if cnt_endyr_APRP>numyrs+1:
                print("error - you are requesting years that do not exist... double check your input")
            ds_cnt[vv]=ds_temp[vv].isel(time=slice((cnt_endyr_APRP-Nyrs_APRP)*12,cnt_endyr_APRP*12)).values
            print(np.shape(ds_cnt[vv])) 
        if Ftype=='default':
            #print(filename_cnt)
            ds_temp=xr.open_dataset(filename_cnt,decode_times=False)
            #print(ds_temp['time'][0],ds_temp['time'][-1])
            ds_temp['time'].values=ds_temp['time'].values-15 
            ds_temp=xr.decode_cf(ds_temp,decode_times=True)
            #print(ds_temp['time'][0],ds_temp['time'][-1])    
            numyrs=int(np.floor(len(ds_temp['time'])/12))
            cnt_endyr_APRP=(cnt_startyr_APRP+Nyrs_APRP-int(cnt_startyr))
            if cnt_endyr_APRP>numyrs+1:
                print("error - you are requesting years that do not exist... double check your input")
            ds_cnt[vv]=ds_temp[vv].isel(time=slice((cnt_endyr_APRP-Nyrs_APRP)*12,cnt_endyr_APRP*12)).values
            print(np.shape(ds_cnt[vv])) 

        ### Load in experiment (ds_exp)
        #print(filename_exp)
        ds_temp=xr.open_dataset(filename_exp,decode_times=False)
        #print(ds_temp['time'][0],ds_temp['time'][-1])
        ds_temp['time'].values=ds_temp['time'].values-15 
        ds_temp=xr.decode_cf(ds_temp,decode_times=True)
        #print(ds_temp['time'][0],ds_temp['time'][-1])    
        numyrs=int(np.floor(len(ds_temp['time'])/12))
        exp_endyr_APRP=(exp_startyr_APRP+Nyrs_APRP-int(exp_startyr))
        if exp_endyr_APRP>numyrs+1:
            print("error - you are requesting years that do not exist... double check your input")
        ds_exp[vv]=ds_temp[vv].isel(time=slice((exp_endyr_APRP-Nyrs_APRP)*12,exp_endyr_APRP*12)).values
        print(np.shape(ds_exp[vv]))    

    lon=ds_temp['lon'].values
    lat=ds_temp['lat'].values
    print("done reading in data")

    return ds_cnt,ds_exp,expnames,lon,lat

def calcgma(swdntoa,swuptoa,swdnsfc,swupsfc):

    Qsd=swdnsfc/swdntoa
    alphaplan=swuptoa/swdntoa
    alphasurf=swupsfc/swdnsfc
    mu = alphaplan + Qsd*(1-alphasurf)
    absorp = 1-mu
    gamma = (mu-Qsd)/(mu-alphasurf*Qsd)

    return gamma,mu,alphasurf

#function to read in all output needed
def calcavg(pert_eachtime,cnt_monthlymean):
    array_to_return=(pert_eachtime+cnt_monthlymean)/2
   # for ii in range(12):
   #     array_to_return[ii::12,:,:]=(exp_eachtime[ii::12,:,:]+cnt_monthlymean[ii,:,:])/2
    
    return array_to_return
        
def docalc(ds_cnt,ds_pert,names):
    #define dictionaries for output
    planetary_albedo={}
    surface_albedo={}
    
    #clear sky, overcast, cloud-sky and other dicts
    gclr={};muclr={};alphasurfclr={}
    goc={};muoc={};alphasurfoc={}
    gcld={};mucld={}
    totcld={}
    
    #defining dicts for outputs
    albedo_clrsky_alpha={};albedo_clrsky_alpha_0={}
    albedo_ocsky_alpha={};albedo_ocsky_alpha_0={}
    albedo_cldsky_mu={};albedo_cldsky_mu_0={}
    albedo_cldsky_gamma={};albedo_cldsky_gamma_0={}
    albedo_cldfrac={};albedo_cldfrac_0={}
    albedo_clrsky_mu={};albedo_clrsky_mu_0={}
    albedo_clrsky_gamma={};albedo_clrsky_gamma_0={}
    
    ds_list={names[0]:ds_cnt, names[1]:ds_pert}
    
    for nn in names:
        ##grab all sky TOA SW
        swdntoa=ds_list[nn]['SOLIN']
        swdnsfc=ds_list[nn]['FSDS']
        swuptoa=swdntoa-ds_list[nn]['FSNT']
        swupsfc=swdnsfc-ds_list[nn]['FSNS']

        ##grab clear sky TOA SW
        swdntoaclr=swdntoa#ds['swdn_toa_clr'].isel(time=sea)
        swuptoaclr=swdntoaclr-ds_list[nn]['FSNTC']
        swdnsfcclr=ds_list[nn]['FSDSC']
        swupsfcclr=swdnsfcclr-ds_list[nn]['FSNSC']
        
        totcld[nn]=ds_list[nn]['CLDTOT'] #CHECK HERE if need to divide by 100%!
        
        #calculate planetary and surface albedo
        planetary_albedo[nn]=swuptoa/swdntoa 
        surface_albedo[nn]=swupsfc/swdnsfc 
        
        #run aprp routine for clear sky
        gclr[nn],muclr[nn],alphasurfclr[nn]=calcgma(swdntoaclr,swuptoaclr,swdnsfcclr,swupsfcclr)

        #calculate overcast sky SW quantities
        swdntoaoc=(swdntoa-(1-totcld[nn])*swdntoaclr)/(totcld[nn])
        swuptoaoc=(swuptoa-(1-totcld[nn])*swuptoaclr)/(totcld[nn])
        swdnsfcoc=(swdnsfc-(1-totcld[nn])*swdnsfcclr)/(totcld[nn])
        swupsfcoc=(swupsfc-(1-totcld[nn])*swupsfcclr)/(totcld[nn])
        
        #run APRP on overcast quantities and calculate cloud-sky
        goc[nn],muoc[nn],alphasurfoc[nn]=calcgma(swdntoaoc,swuptoaoc,swdnsfcoc,swupsfcoc) 
        mucld[nn]=muoc[nn]/muclr[nn]
        gcld[nn]=1-(1-goc[nn])/(1-gclr[nn])
        
    #using method 12a to calculate differences - might be nice to implement 12b or 12d later
    cn=names[0]
    pn=names[1]
    #for nn in names[1:]:
    mumeanclr=calcavg(muclr[pn],muclr[cn])
    gammameanclr=calcavg(gclr[pn],gclr[cn])
    totmean=calcavg(totcld[pn],totcld[cn])
    mumeanoc=calcavg(muoc[pn],muoc[cn])
    gammameanoc=calcavg(goc[pn],goc[cn])
    gmeancld=calcavg(gcld[pn],gcld[cn])
    alphameanoc=calcavg(alphasurfoc[pn],alphasurfoc[cn])
    mumeancld=calcavg(mucld[pn],mucld[cn])
    alphameanclr=calcavg(alphasurfclr[pn],alphasurfclr[cn])
    gammameancld=calcavg(gcld[pn],gcld[cn])
        
    #Albedo from clearsky alpha 
    Aclr=mumeanclr*alphasurfclr[pn]*(1-gammameanclr)**2/(1-alphasurfclr[pn]*gammameanclr)
    Aclr0=mumeanclr*alphasurfclr[cn]*(1-gammameanclr)**2/(1-alphasurfclr[cn]*gammameanclr)
    albedo_clrsky_alpha=(Aclr)*(1-totmean)  
    albedo_clrsky_alpha_0=Aclr0*(1-totmean)
        
    #Albedo from overcast alpha
    Aoc=mumeanoc*alphasurfoc[pn]*(1-gammameanoc)**2/(1-alphasurfoc[pn]*gammameanoc)
    Aoc0=mumeanoc*alphasurfoc[cn]*(1-gammameanoc)**2/(1-alphasurfoc[cn]*gammameanoc)
    albedo_ocsky_alpha=(Aoc)*totmean
    albedo_ocsky_alpha_0=Aoc0*totmean   
        
    #Albedo from cloudy sky mu
    part1muc=mumeanclr*mucld[pn]*gammameanoc
    part2muc=alphameanoc*mumeanclr*mucld[pn]*(1-gammameanoc)**2/(1-alphameanoc*gammameanoc)            
    part1muc0=mumeanclr*mucld[cn]*gammameanoc
    part2muc0=alphameanoc*mumeanclr*mucld[cn]*(1-gammameanoc)**2/(1-alphameanoc*gammameanoc)            
    albedo_cldsky_mu=((part1muc+part2muc))*totmean
    albedo_cldsky_mu_0=(part1muc0+part2muc0)*totmean
        
    #Albedo from cloudy sky gamma
    part1gc=mumeanoc*(1-(1-gammameanclr)*(1-gcld[pn]))  
    part2gc=alphameanoc*mumeanoc*((1-gammameanclr)*(1-gcld[pn]))**2/(1-alphameanoc*(1-(1-gammameanclr)*(1-gcld[pn])))
    part1gc0=mumeanoc*(1-(1-gammameanclr)*(1-gcld[cn]))  
    part2gc0=alphameanoc*mumeanoc*((1-gammameanclr)*(1-gcld[cn]))**2/(1-alphameanoc*(1-(1-gammameanclr)*(1-gcld[cn])))
    albedo_cldsky_gamma= (part1gc+part2gc)*totmean
    albedo_cldsky_gamma_0= (part1gc0+part2gc0)*totmean
        
    #Albedo from cloud amount      
    part1cc=(1-totcld[pn])*(mumeanclr*gammameanclr+alphameanclr*mumeanclr*(1-gammameanclr)**2/(1-alphameanclr*gammameanclr))
    part2cc=(totcld[pn])*(mumeanoc*gammameanoc+alphameanoc*mumeanoc*(1-gammameanoc)**2/(1-alphameanoc*gammameanoc))
    part1cc0=(1-totcld[cn])*(mumeanclr*gammameanclr+alphameanclr*mumeanclr*(1-gammameanclr)**2/(1-alphameanclr*gammameanclr))
    part2cc0=(totcld[cn])*(mumeanoc*gammameanoc+alphameanoc*mumeanoc*(1-gammameanoc)**2/(1-alphameanoc*gammameanoc))
    albedo_cldfrac=part1cc+part2cc
    albedo_cldfrac_0=part1cc0+part2cc0
        
    #Albedo from clear sky mu
    part1muc = (1-totmean)*(muclr[pn]*gammameanclr+alphameanclr*muclr[nn]*(1-gammameanclr)**2/(1-alphameanclr*gammameanclr))
    part2muc = totmean*(muclr[pn]*mumeancld*gammameanoc + alphameanoc*muclr[pn]*mumeancld*(1-gammameanoc)**2/(1-alphameanoc*gammameanoc))
    part1muc0 = (1-totmean)*(muclr[cn]*gammameanclr+alphameanclr*muclr[cn]*(1-gammameanclr)**2/(1-alphameanclr*gammameanclr))
    part2muc0 = totmean*(muclr[cn]*mumeancld*gammameanoc + alphameanoc*muclr[cn]*mumeancld*(1-gammameanoc)**2/(1-alphameanoc*gammameanoc))    
    albedo_clrsky_mu=part1muc+part2muc 
    albedo_clrsky_mu_0= part1muc0+part2muc0
        
    #Albedo from clear sky gamma
    part1gcl = (1-totmean)* (mumeanclr*gclr[pn] +alphameanclr*mumeanclr*(1-gclr[pn])**2 /(1-alphameanclr*gclr[pn])) 
    part2gcl = totmean*(mumeanoc*(1-(1-gclr[pn])*(1-gammameancld)) )
    part3gcl = totmean*(alphameanoc*mumeanoc*((1-gclr[pn])*(1-gammameancld))**2/(1-alphameanoc*(1-(1-gclr[pn])*(1-gammameancld))))           
    part1gcl0 = (1-totmean)* (mumeanclr*gclr[cn] +alphameanclr*mumeanclr*(1-gclr[cn])**2 /(1-alphameanclr*gclr[cn])) 
    part2gcl0 = totmean*(mumeanoc*(1-(1-gclr[cn])*(1-gammameancld)) )
    part3gcl0 = totmean*(alphameanoc*mumeanoc*((1-gclr[cn])*(1-gammameancld))**2/(1-alphameanoc*(1-(1-gclr[cn])*(1-gammameancld))) )           
    albedo_clrsky_gamma=(part1gcl+part2gcl+part3gcl)
    albedo_clrsky_gamma_0=part1gcl0+part2gcl0+part3gcl0
        
    #one big dictionary to make easier
    returnme={}
    returnme['planetary_albedo']=planetary_albedo
    returnme['surface_albedo']=surface_albedo
    returnme['albedo_clrsky_alpha']=albedo_clrsky_alpha
    returnme['albedo_clrsky_alpha_0']=albedo_clrsky_alpha_0
    returnme['delta_albedo_clrsky_alpha']=albedo_clrsky_alpha-albedo_clrsky_alpha_0
    returnme['albedo_ocsky_alpha']=albedo_ocsky_alpha
    returnme['albedo_ocsky_alpha_0']=albedo_ocsky_alpha_0
    returnme['delta_albedo_ocsky_alpha']=albedo_ocsky_alpha-albedo_ocsky_alpha_0
    returnme['albedo_cldsky_mu']=albedo_cldsky_mu
    returnme['albedo_cldsky_mu_0']=albedo_cldsky_mu_0
    returnme['delta_albedo_cldsky_mu']=albedo_cldsky_mu-albedo_cldsky_mu_0
    returnme['albedo_cldsky_gamma']=albedo_cldsky_gamma
    returnme['albedo_cldsky_gamma_0']=albedo_cldsky_gamma_0
    returnme['delta_albedo_cldsky_gamma']=albedo_cldsky_gamma-albedo_cldsky_gamma_0
    returnme['albedo_cldfrac']=albedo_cldfrac
    returnme['albedo_cldfrac_0']=albedo_cldfrac_0
    returnme['delta_albedo_cldfrac']=albedo_cldfrac-albedo_cldfrac_0
    returnme['albedo_clrsky_mu']=albedo_clrsky_mu
    returnme['albedo_clrsky_mu_0']=albedo_clrsky_mu_0
    returnme['delta_albedo_clrsky_mu']=albedo_clrsky_mu-albedo_clrsky_mu_0
    returnme['albedo_clrsky_gamma']=albedo_clrsky_gamma
    returnme['albedo_clrsky_gamma_0']=albedo_clrsky_gamma_0
    returnme['delta_albedo_clrsky_gamma']=albedo_clrsky_gamma-albedo_clrsky_gamma_0
    return returnme

