# -*- coding: utf-8 -*-
__author__='minqiang'
# this file is created for read MIPAS data for CFC validation
# 

import logging,ftir
from ftir import getlogger
rootlogger=logging.getLogger(__name__)
from numpy import *
import coda,h5py,pygrib,os,subprocess,sys,shutil,pickle,glob,urllib2
from datetime import datetime,timedelta
try: from pyhdf.SD import SD,SDC
except: pass
from ftir.tools import hdf_store,hdf5todict
import matplotlib.pyplot as plt
import os
import ftir.geoms as g
from netCDF4 import Dataset
from tools_min import D_for_boundext


def get_colocated_data(mol='',lat=-21,lon=55,co_locate='5x5',min_alt=7.0,max_alt=20.0,hdf=False,outpath='/home/minqiang/working/sfit4/test/mipas/',versions=None,logger=rootlogger):
  """get the co-located MIPAS data for a specific latitude and logitude  and also calculate the partical column for some comparison
  mol : <str> the target mol ,such as F-11
  lat : <float> latitude
  lon : <float> longitude
  min_alt : <float> the bottom of partial column
  max_alt : <float> the top of partial column
  hdf : bool  store into a hdf file or not
  co_locate: the co-located mode you choose
  Note that: some products have different versions! versions=None,['V3O','V5H','V5R']...  """
  logger=getlogger(logger,'get_colocated_data')  
  if co_locate=='10x30': logger.info('the co-located box is chooese as 10 x 30'); dlat=10;dlon=30  
  if co_locate=='5x15': logger.info('the co-located box is chooese as 5 x 15'); dlat=5;dlon=15
  if co_locate=='5x5': logger.info('the co-located box is chooese as 5 x 5'); dlat=5;dlon=5   
  if co_locate=='2x2': logger.info('the co-located box is chooese as 2 x 2'); dlat=2;dlon=2 
  if co_locate=='2x5': logger.info('the co-located box is chooese as 2 x 5'); dlat=2;dlon=5 
  
  if not mol : logger.error('No mol is provide!'); return 1
  if mol=='CCl3F': mol='F-11'
  if mol=='CCl2F2': mol='F-12'
  if mol=='CHF2Cl': mol='F-22'
  mipasfpath='/data/SATELLITE/ENVISAT/MIPAS/L2/KIT/archive/share.lsdf.kit.edu/imk/asf/sat/mipas-export/Data_by_Target' # this is the KIT version <if you want to select other products, change this folder>
  mipasfiles=[]
  for e in os.walk(os.path.join(mipasfpath,mol)):
    if len(e[1])==0 and len(e[2])>0:
      if not versions: # choose the latest version
	mipasfiles.append(os.path.join(e[0],e[2][-1]))
      else:  # only selected the specify version
	if type(versions)==str:versions=[versions]
	for efile in e[2]:
	  for version in versions:
	    if version in efile:mipasfiles.append(os.path.join(e[0],efile))
  if not mipasfiles : logger.error('There is no files exits in %s'%os.path.join(mipasfpath,mol.upper()))
  out={};out['lat']=[];out['lon']=[];out['target']=[];out['pressure']=[];out['temperature']=[]
  out['altitude']=[];out['target_noise']=[];out['time']=[];out['target_pc']=[];out['target_pc_error']=[]
  out['target_tc']=[];out['target_tc_error']=[];out['avk_dia']=[]
  j2000_jd=timedelta(2451545);j1=datetime(2000,1,1,12)
  for i,mipasfile in enumerate(mipasfiles):
    mipas_nc_file=mipasfile
    fh=Dataset(mipas_nc_file,mode='r')
    latitude=fh.variables['latitude']
    longitude=fh.variables['longitude']
    logger.info('read file : %s'%mipasfile)
    for j in arange(len(latitude)):
      if abs(latitude[j]-lat) > dlat : continue
      if abs(longitude[j]-lon) > dlon : continue
      out['lat'].append(latitude[j])
      out['lon'].append(longitude[j])
      if fh.variables['time'].units == 'julian days':out['time'].append(j1+timedelta(fh.variables['time'][j])-j2000_jd)	  
      elif 'days since' in fh.variables['time'].units : year,month,day=fh.variables['time'].units.split()[2].split('-');out['time'].append(datetime(int(year),int(month),int(day))+timedelta(fh.variables['time'][j]))
      else: logger.error('error with the time struct')
      out['target'].append(fh.variables['target'][:,j])
      out['target_noise'].append(fh.variables['target_noise_error'][:,j])
      out['pressure'].append(fh.variables['pressure'][:,j])
      out['temperature'].append(fh.variables['temperature'][:,j])
      out['altitude'].append(fh.variables['altitude'][:,j])
      out['avk_dia'].append(fh.variables['akm_diagonal'][:,j])
      alt=fh.variables['altitude'][:,j];tem=fh.variables['temperature'][:,j];pre=fh.variables['pressure'][:,j]
      pc=[];pc_error=[];boundary_top=[];boundary_bottom=[]
      for n in arange(len(alt)):
	if n==0: air_dry=pre[n]*1000*(alt[n+1]-alt[n])/2.0/8.314/tem[n];boundary_top.append((alt[n+1]+alt[n])/2.0);boundary_bottom.append(alt[n])
	elif n==len(alt)-1: air_dry=pre[n]*1000*(alt[n]-alt[n-1])/2.0/8.314/tem[n];boundary_top.append(alt[n]);boundary_bottom.append((alt[n]+alt[n-1])/2.0)
	else: air_dry=pre[n]*1000*(alt[n+1]-alt[n-1])/2.0/8.314/tem[n];boundary_top.append((alt[n+1]+alt[n])/2.0);boundary_bottom.append((alt[n]+alt[n-1])/2.0)
	pc.append(air_dry*fh.variables['target'][n,j])
	pc_error.append(air_dry*fh.variables['target_noise_error'][n,j])
      boundary=[boundary_bottom,boundary_top]
      weighting_for_pc=D_for_boundext(array(boundary),array([[min_alt],[max_alt]]))
      weighting_for_tc=D_for_boundext(array(boundary),array([[0.0],[120.0]]))
      #print weighting_for_pc,weighting_for_tc
      #print array(boundary)
      gas_patical_column=sum(array(pc)*weighting_for_pc)   
      gas_patical_column_error=sum(array(pc_error)*weighting_for_pc)	
      gas_total_column=sum(array(pc)*weighting_for_tc)   
      gas_total_column_error=sum(array(pc_error)*weighting_for_tc)
      out['target_pc'].append(gas_patical_column)
      out['target_pc_error'].append(gas_patical_column_error)   
      out['target_tc'].append(gas_total_column)
      out['target_tc_error'].append(gas_total_column_error)   
  if not out : logger.error('there is no co-located MIPAS data for %s'%mol); return 1
  if hdf : 
    for key in out.keys():
      out[key]=array(out[key])  # transfer to the array
    with h5py.File(outpath+mol+'.hdf','w') as hdfid: ftir.tools.hdf_store(hdfid,'',out) 
  return out;
  
