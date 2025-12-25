#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 15:25:03 2025

@author: joris
"""
import glob
import numpy as np
#import random
import pandas as pd
#import datetime
from datetime import date,time, timedelta, datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#from oct2py import octave
import scipy.ndimage
import sys
import os
import sys
import re
from datetime import date,time, timedelta

version='dev version'
red='\x1b[41m'
orange='\x1b[43m'
green='\x1b[42m'
black='\x1b[40m'
black='\x1b[47m'

print('Loading GlyApp '+version)
# Get the encoding of the file in relation to the file type
FILE_ENCODING = {
	'capteur_medtronics': 'utf-16',
	'capteur_ipro': 'utf-16le',
	'capteur_standard': 'utf-8',
	'capteur_freestyle_pro': 'utf-8',
	'capteur_freestyle': 'utf-8',
	'capteur_dexcom_canada_mmol': 'utf-8',
	'capteur_medtronics_canada_mmol': 'utf-8',
	'accelero_gt3x': 'ASCII',
	'accelero_gt1m': 'ASCII',
	'alimentation': 'iso8859',
	'insulinelente': 'iso8859',
	'insulinerapide': 'iso8859',
	'insuline': 'utf-8',
	'sport': 'iso8859',
}

def test_all():
	# test all folder in DATA
	base='./Data/'
	folders=glob.glob(base+'*')
	score=0
	for f in folders:
		try:
			read_Glu(f[len(base):])
			score+=1
		except:
			print(red+'Error reading patient'+f +'\n'+black)
	print('\n \n ==== \n SCORE : ',score/len(folders)*100,' % of successful reading !')
	return None
	
def read_libraries(delimiter=',',encoding='utf-8'):
	sensors=pd.read_csv('sensor_library.csv',delimiter=delimiter,encoding=encoding)
	filenames=pd.read_csv('filename_sensors.csv',delimiter=delimiter)
	return filenames,sensors

def import_libraries_pyodide():
	url = 'https://docs.google.com/spreadsheets/d/1KH61giiVdZmQJ8h97awgeSCJKCYvDwdM6-nsSYvjSxc/gviz/tq?tqx=out:csv&sheet={sensor_library}'
	with open('sensor_library.csv','w') as f:
		f.write(pyodide.http.open_url(url).read())
	url = 'https://docs.google.com/spreadsheets/d/1Ya7Y9joet36BuhjC7pecI9VnejrPgPFFIm6Z4mo29IU/gviz/tq?tqx=out:csv&sheet={sensor_library}'
	with open('filename_sensors.csv','w') as f:
		f.write(pyodide.http.open_url(url).read())
	return None

def conversion_factor(units):
	factor=1

def import_libraries_urllib():
	import urllib
	url = 'https://docs.google.com/spreadsheets/d/1KH61giiVdZmQJ8h97awgeSCJKCYvDwdM6-nsSYvjSxc/gviz/tq?tqx=out:csv&sheet={sensor_library}'
	urllib.request.urlretrieve(url,'sensor_library.csv')
	url = 'https://docs.google.com/spreadsheets/d/1Ya7Y9joet36BuhjC7pecI9VnejrPgPFFIm6Z4mo29IU/gviz/tq?tqx=out:csv&sheet={sensor_library}'
	urllib.request.urlretrieve(url,'filename_sensors.csv')
	return None

def conversion_factor(units):
	factor=1
	if units=='mmol/l':
		factor= 18.02  # conversion factor mmol/L to mg/L
	return factor

def date_col_num(d):
	df = pd.DataFrame(d, columns=['year','month','day','hour','minute','second'])
	num=datenum(df)
	return num

def sort_times(GluTime,GluValue):
	# ===== Sort by ascending time =====
	GluTime = np.array(GluTime)
	GluValue = np.array(GluValue)
	
	# Possibly reorder year month day
	if GluTime.shape[0]>0:
		GluTime_temp=np.copy(GluTime)
		idyear=np.where(GluTime[0,:]>1000)[0][0]
		GluTime_temp[:,0]=GluTime[:,idyear]
		GluTime_temp[:,idyear]=GluTime[:,0]
		GluTime=GluTime_temp
	
	datetimes = [datetime(*map(int, t)) for t in GluTime]
	sorted_idx = np.argsort(datetimes)
	GluTime = GluTime[sorted_idx]
	GluValue = GluValue[sorted_idx]
	return GluTime,GluValue


def clean_data(GluTime,GluValue):
	# ===== remove nan data =====
	GluTime = np.array(GluTime)
	GluValue = np.array(GluValue)
	# remove nan and years that are = lower than 2000
	idnan=np.where(np.isnan(GluValue[:,0])|(GluTime[:,0]==2000))[0]
	return np.delete(GluTime,idnan,axis=0),np.delete(GluValue,idnan,axis=0)

def datenum(d):
	dates = pd.to_datetime(d)
	#print(type(dates))
	#res=pd.to_numeric((dates - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s'))/3600/24.
	if isinstance(dates, pd.Timestamp):
		res1=((dates - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s'))/3600/24.
		#print(type(dates),res1)
	else:
		if len(dates)>1: # to avoid warning
			res1=np.float64((dates - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s'))/3600/24.
		else:
			res1=np.float64((dates.iloc[0] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s'))/3600/24.
		#print(type(dates),res1)
	#print(d,'res:',res.to_numpy(),'res1:',res1)
	#print(res1)
	return res1

def numdate(n):
	return pd.to_datetime((n*3600*24)*pd.Timedelta('1s')+pd.Timestamp("1970-01-01"))


def calc_MAGE(Glu_interp,Time,nstd,Tav):
	# MAGE
	Glu_std=np.nanstd(Glu_interp)
	
	#Glu_f=Gaussian1DKernel(Glu,Tav) # Smoothing Glucose over Tav minutes
	Glu_f=scipy.ndimage.gaussian_filter1d(Glu_interp, 5)
	Glu_max=scipy.ndimage.maximum_filter(Glu_f, size=int(Tav))
	idmax=np.where(Glu_max==Glu_f)[0]
	Glu_min=scipy.ndimage.minimum_filter(Glu_f, size=int(Tav))
	idmin=np.where(Glu_min==Glu_f)[0]
	#plt.plot(Glu_max);plt.plot(Glu_f);plt.plot(idmax,Glu_interp[idmax],'ko');plt.plot(idmin,Glu_interp[idmin],'ro')
	
	# Find local Peaks
	
	idAll=np.sort(np.hstack((idmax,idmin)));
	cutoff=nstd; #(in SD)
	idAll_filt=idAll;
	
	# Compute Peak amplitude 
	I=1;
	while (I==1)&(len(idAll_filt)>=2):
		#print(len(idAll_filt))
		Amplitude=(np.abs(Glu_interp[idAll_filt[2:]]-Glu_interp[idAll_filt[1:-1]])+np.abs(Glu_interp[idAll_filt[1:-1]]-Glu_interp[idAll_filt[:-2]]))/2;
		Amplitude=np.hstack((np.abs(Glu_interp[idAll_filt[1]]-Glu_interp[idAll_filt[0]]),Amplitude,np.abs(Glu_interp[idAll_filt[-1]]-Glu_interp[idAll_filt[-2]])));
		# Sort amplitudes from smallest to largest
		idsort=np.argsort(Amplitude)
		A_sort=Amplitude[idsort]
		# 2/ Remove smallest amplitude under cutoff*STD_glu
		#print(A_sort)
		idcut=np.where(A_sort<cutoff*Glu_std)[0]
		if len(idcut)>0:
			idAll_filt=np.delete(idAll_filt,idsort[idcut[0]]);
		else:
			I=0
	#plt.plot(Glu_max);plt.plot(Glu_f);plt.plot(idAll_filt,Glu_interp[idAll_filt],'ko')
	# 3/ Remove Non-switching extremas
	idbad=np.where((Glu_interp[idAll_filt[2:]]-Glu_interp[idAll_filt[1:-1]])*
														(Glu_interp[idAll_filt[1:-1]]-Glu_interp[idAll_filt[:-2]])>0)[0]
	idAll_filt=np.delete(idAll_filt,idbad+1)
	# 4/ Difference between peaks
	Amplitude_real=(np.diff(Glu_interp[idAll_filt]));
	# Average time between two peaks
	Time_Amplitude=(Time[idAll_filt[1:-1]]+Time[idAll_filt[2:]])/2;
	return Amplitude_real,Time_Amplitude,idAll_filt


def CONGAn(Glu_interp):
	# # CONGAn, computed from start_Hour to start_Hour each day
	CONGAn=[]
	# if nDays>0:
	# 	for n in range(1,15): # Hours of shift
	# 		ndt=60*n; #Equivalent in time resolution
	# 		Glu_temp=np.hstack((Glu_interp[ndt:],np.zeros(ndt)+np.nan));
	# 		Dif=Glu_temp-Glu_interp;
	# 		CONGAn_all=Dif.reshape(24*60,nDays);
	# 		CONGAn_all[-ndt:,:]=np.nan; # delete the last ndt daily values that concerns the following day
	# 		CONGAn.append(np.nanstd(CONGAn_all,axis=0))
	# 	CONGA_n_mean=np.nanmean(np.array(CONGAn),axis=1);
	# 	CONGA_n_mean_err=np.nanstd(np.array(CONGAn),axis=1)/np.sqrt(nDays)
	# else:
	# 	CONGA_n_mean=np.nan
	# 	CONGA_n_mean_err=np.nan
	#if nDays>0:
	for n in range(1,9): # Hours of 
		try:
			ndt=60*n; #Equivalent in time resolution
			Glu_temp=np.hstack((Glu_interp[ndt:],np.zeros(ndt)+np.nan));
			Dif=Glu_temp-Glu_interp;
			CONGAn.append(np.nanstd(Dif))
		except:
			CONGAn.append(np.nan)
	return np.array(CONGAn)

# Calculation of hypoglycemias
def calc_hypo(Glu_interp,T_interp,time_hypo,glu_hypo):
# =============================================================================
# 	# Normal hypoglycemia
# 	time_ep_start=15; # Min duration of Hypo episod  in seconds(start)
# 	time_ep_stop=15; # Min duration of Hypo episod  (end)
# 	hypo_start=54;
# 	hypo_stop=70;
# 	# Prolongated hypoglycemia
# 	time_ep_start_pro=120; # Min duration of Hypo episod  in seconds(start)
# 	time_ep_stop_pro=120; # Min duration of Hypo episod  (end)
# 	hypo_start_pro=54;
# 	hypo_stop_pro=54;
# =============================================================================
	hypo_start,hypo_stop=glu_hypo
	time_ep_start,time_ep_stop=time_hypo
	Nreading_start=time_ep_start;
	Nreading_stop=time_ep_stop;
	Start_hypo=[];
	Stop_hypo=[];
	
	# normal hypo
	Nhypo=0; # nb of hyporheic ep
	i=1;
	while i<(Glu_interp.shape[0]-Nreading_start):
		#check if hypoglycemic start
		start=np.all(Glu_interp[i:i+Nreading_stop]<hypo_start);
		if start==1: # hypoglycemic start
			Nhypo=Nhypo+1;
			Start_hypo.append(T_interp[i]);
			stop=1;
			i=i+Nreading_start+1;
			while (stop) & (i<(Glu_interp.shape[0]-Nreading_stop)):
				stop=np.any(Glu_interp[i:i+Nreading_stop]<hypo_stop);
				i=i+1;
			Stop_hypo.append(T_interp[i-1]);
		else:
			i=i+1;
	
	Start_hypo=np.array(Start_hypo)
	Stop_hypo=np.array(Stop_hypo)
	return Start_hypo,Stop_hypo
# ===============INITIALISATION================================================



def plot_patient(patient='GZ2',encoding='utf-8',
				 start=None,
				 end=None,
				 start_Hour=7,
				 hypo=70,hyper=180,
				 superhypo=54,superhyper=250,
				 plot_acc=False,
				 filt_acc=1,
				 seuils_acc=[0,100,250,500,1000,2000],
 				 plot_temp=0,
				 plot_bpm=False,
				 filt_bpm=1,
				 seuils_bpm=[20,300],
				 plot_alt=0,
				 savefig=True,
				 lw=2):
	
	xformatter = mdates.DateFormatter('%H:%M')
	GluTime_all,GluValue_all=read_Glu(patient,encoding=encoding)
	T_interp_all=np.arange(date_col_num(GluTime_all[0,:].reshape(1,-1)),date_col_num(GluTime_all[-1,:].reshape(1,-1)),1/(60*24)) # reinterpolate data on 1 min intervals
	Glu_interp_all=np.interp(T_interp_all,date_col_num(GluTime_all),GluValue_all[:,1])
	Glu_all=np.copy(GluValue_all)
	# Date as a panda format
	GluTime_all_pd=numdate(date_col_num(GluTime_all))
	#Index_all_pd=np.arange(len(GluTime_pd))
	
	# Read accelero :
	if plot_acc:
		try:
			accelero=read_accelero(patient)
			# filter
			#accelero['Norm']=scipy.ndimage.median_filter(accelero['Norm'].to_numpy(),nfilt)
			#N=scipy.ndimage.gaussian_filter1d(accelero['Norm'].to_numpy(),nfilt)
			Nini=scipy.ndimage.median_filter(accelero['Norm'].to_numpy(),
			filt_acc)
			N=np.copy(Nini)
			# si seuils, modifier N
			print(N.max(),N.min())
			for i in range(len(seuils_acc)-1):
				N[(Nini>seuils_acc[i])&(Nini<seuils_acc[i+1])]=i
		except:
			print('Warning: no accelerometer was read')
			plot_acc=False

	# Read cardio :
	if (plot_bpm)|(plot_alt>0)|(plot_temp>0):
		try:
			cardio=read_cardio(patient)
			# filter
			nfilt=50
			#accelero['Norm']=scipy.ndimage.median_filter(accelero['Norm'].to_numpy(),nfilt)
			#N=scipy.ndimage.gaussian_filter1d(accelero['Norm'].to_numpy(),nfilt)
			BPMini=scipy.ndimage.median_filter(cardio['HR (bpm)'].to_numpy(),filt_bpm)
			BPM=np.copy(BPMini)
			#print(np.nanmin(BPM),np.nanmax(BPM))
			for i in range(len(seuils_bpm)-1):
				#print(seuils_bpm[i])
				BPM[(BPMini>seuils_bpm[i])&(BPMini<=seuils_bpm[i+1])]=i
			#print(np.nanmin(BPM),np.nanmax(BPM))
			ALT=scipy.ndimage.median_filter(cardio['Altitude (m)'].to_numpy(),nfilt)*plot_alt
			TMP=(scipy.ndimage.median_filter(cardio['Temperatures (C)'].to_numpy(),nfilt)+5)*plot_temp
		except:
			print('Warning: no cardio or temp was read')
			plot_bpm=False
			plot_alt=0
			plot_temp=0

	# Intervals
	if not start:
		start = numdate(date_col_num(GluTime_all[:1,:]))
	else:
		start = pd.Timestamp(start)
	if not end:
		end = numdate(date_col_num(GluTime_all[-1:,:]))
	else:
		end = pd.Timestamp(end)
	
	if start.hour<7:
		start=start+pd.Timedelta(days=-1)
	start_1=start.replace(hour=start_Hour,minute=0,second=0)
	if end.hour>7:
		end=end+pd.Timedelta(days=1)
	end_1=end.replace(hour=start_Hour,minute=0,second=0)
	
	print(start,end)
	ndays=(end_1-start_1).days
	ndays=np.maximum(np.minimum(ndays,7),2) # maximum 7 days is plotted
	
	fig,ax=plt.subplots(ndays,1,figsize=(15,ndays*0.2*15),sharex=True)
	for n in range(ndays):
# =============================================================================
		# Plot Accelero
# =============================================================================
		if plot_acc:
			isin_pd=(accelero['Date']>start_1+pd.Timedelta(days=n))&(accelero['Date']<start_1+pd.Timedelta(days=n+1))
			img=N[isin_pd].reshape(1,-1)
			#print(n,len(isin_pd))
			ax[n].imshow(img,cmap=plt.cm.YlOrBr,extent=[start_1,start_1+pd.Timedelta(days=1),0,superhyper*1.2],alpha=0.4,aspect='auto')

		if plot_bpm:
			isin_pd=(cardio['Time']>start_1+pd.Timedelta(days=n))&(cardio['Time']<start_1+pd.Timedelta(days=n+1))
# 			ny=200
# 			img=BPM[isin_pd].reshape(1,-1)
# 			if len(seuils_bpm)==0:
# 				img=BPM[isin_pd].reshape(-1,1)
# 				img=np.tile(img,ny).T
# 				x,y=np.meshgrid(np.arange(img.shape[1]),np.arange(img.shape[0]))
# 				std=ny/800*(BPM[isin_pd]-30).reshape(1,-1)
# 				img=img*np.exp(-(y-ny/2)**2/(2*std**2))
			#print(BPM[isin_pd])
			Timetemp=cardio['Time'][isin_pd].copy().reset_index()
			ax[n].fill_between(Timetemp['Time']-pd.Timedelta(days=n),150-(BPM[isin_pd]-50)*1,150+(BPM[isin_pd]-50)*1,color='Red',alpha=0.5,edgecolor='w')
# 			for i in range(len(seuils_bpm)-1):
# #				print(BPM[isin_pd],len(BPM[isin_pd]))
# 				id_seuil=np.where((BPM[isin_pd]>seuils_bpm[i])&((BPM[isin_pd]<=seuils_bpm[i+1])))[0]
# 				print(len(id_seuil))
# # 				print(cardio['Time'][isin_pd])
# # 				print(id_seuil,cardio['Time'][isin_pd])
# 				ax[n].fill_between(Timetemp['Time'][id_seuil]-pd.Timedelta(days=n),150-(BPM[isin_pd][id_seuil]-50)*1,150+(BPM[isin_pd][id_seuil]-50)*1,color=plt.cm.PuRd((i+1)/(len(seuils_bpm)-1.)),alpha=0.5,edgecolor='w')
			#print(BPM[isin_pd])
			#print(n,len(isin_pd))
			#ax[n].imshow(img,cmap=plt.cm.PuRd,extent=[start_1,start_1+pd.Timedelta(days=1),0,superhyper*1.2],alpha=0.4,aspect='auto')

		if plot_alt>0:
			isin_pd=(cardio['Time']>start_1+pd.Timedelta(days=n))&(cardio['Time']<start_1+pd.Timedelta(days=n+1))
			#print(ALT[isin_pd])
			ax[n].fill_between(cardio['Time'][isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd)),ALT[isin_pd],color='Gray',alpha=0.4,edgecolor='w')
			alt2gly=lambda x:x/plot_alt
			gly2alt=lambda x:x*plot_alt
			secax = ax[n].secondary_yaxis('right', functions=(alt2gly, gly2alt),color='Gray')
			secax.set_xlabel('Altitude [m]')
		#ax[n].xaxis.grid(True, which='minor')
		
		if plot_temp>0:
			isin_pd=(cardio['Time']>start_1+pd.Timedelta(days=n))&(cardio['Time']<start_1+pd.Timedelta(days=n+1))
			#print(TMP[isin_pd])
			ax[n].plot(cardio['Time'][isin_pd]-pd.Timedelta(days=n),TMP[isin_pd],color='Red',alpha=1)
			tmp2gly=lambda x:x/plot_temp-5
			gly2tmp=lambda x:(x+5)*plot_temp
			secax = ax[n].secondary_yaxis('right', functions=(tmp2gly, gly2tmp),color='Red')
			secax.set_xlabel('Temperature [C]',color='Red')
		#ax[n].xaxis.grid(True, which='minor')
# =============================================================================
# 		Plot Glycemia
# =============================================================================
		# hypo/hyper as areas
		isin_pd=(GluTime_all_pd>start_1+pd.Timedelta(days=n))&(GluTime_all_pd<start_1+pd.Timedelta(days=n+1))
# 		#ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hypo,np.zeros(sum(isin_pd))+hyper,color='g',alpha=0.1)
# 		# Hyper glycemia
# 		glu_temp=np.copy(GluValue_all[isin_pd,1])
# 		#print(np.sum(glu_temp<hyper),glu_temp)
# 		glu_temp[glu_temp<hyper]=np.nan
# 		ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hyper,glu_temp,color='r',alpha=0.5)
# 		# Hypo glycemia
# 		glu_temp=np.copy(GluValue_all[isin_pd,1])
# 		glu_temp[glu_temp>hypo]=np.nan
		T=[start_1,start_1+pd.Timedelta(days=1)]
# 		ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hypo,glu_temp,color='b',alpha=0.5)
		ax[n].plot(T,np.zeros(len(T))+superhyper,'--',linewidth=lw,color='crimson')
		ax[n].plot(T,np.zeros(len(T))+superhypo,'--',linewidth=lw,color='crimson')
		ax[n].plot(T,np.zeros(len(T))+hyper,'-.',color='royalblue',linewidth=lw)
		ax[n].plot(T,np.zeros(len(T))+hypo,'-.',color='royalblue',linewidth=lw)
		ax[n].plot(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),GluValue_all[isin_pd,1],'-',linewidth=lw,color='royalblue')
		ax[n].xaxis.set_major_formatter(xformatter)
		ax[n].set_ylim([GluValue_all[:,1].min(),GluValue_all[:,1].max()])
		ax[n].set_ylim([0,superhyper*1.2])
		ax[n].grid(visible=True)
		ax[n].set_ylabel('Glycemia [mg/dL]',color='royalblue')
		ax[n].tick_params(labelcolor='royalblue')
		ax[n].text(start_1,superhyper*1.25,str(GluTime_all_pd[isin_pd][0].strftime('%A %d-%m-%Y')))
	ax[0].set_xlim([start_1,start_1+pd.Timedelta(days=1)])
	
	ax[n].set_xlabel('Time')
	if savefig:
		save_folder='./Patients/'+patient+'/'
		if not(os.path.isdir(save_folder)):
			os.mkdir(save_folder)
		plt.savefig(save_folder+'Results_'+patient+'.pdf',bbox_inches='tight')

def short_plot_patient(patient='GZ2',encoding='utf-8',
				 start=None,
				 end=None,
				 hypo=70,hyper=180,
				 superhypo=54,superhyper=250,
				 plot_acc=False,
				 filt_acc=1,
				 seuils_acc=[0,100,250,500,1000,2000],
				 plot_bpm=False,
				 filt_bpm=1,
				 seuils_bpm=[20,60,100,120,150,200,300],
				 plot_alt=0,
				 savefig=True,
				 lw=2):

	xformatter = mdates.DateFormatter('%H:%M')
	GluTime_all,GluValue_all=read_Glu(patient,encoding=encoding)
	T_interp_all=np.arange(date_col_num(GluTime_all[0,:].reshape(1,-1)),date_col_num(GluTime_all[-1,:].reshape(1,-1)),1/(60*24)) # reinterpolate data on 1 min intervals
	Glu_interp_all=np.interp(T_interp_all,date_col_num(GluTime_all),GluValue_all[:,1])
	Glu_all=np.copy(GluValue_all)
	# Date as a panda format
	GluTime_all_pd=numdate(date_col_num(GluTime_all))
	#Index_all_pd=np.arange(len(GluTime_pd))
	
	# Read accelero :
	if plot_acc:
		try:
			accelero=read_accelero(patient)
			# filter
			#accelero['Norm']=scipy.ndimage.median_filter(accelero['Norm'].to_numpy(),nfilt)
			#N=scipy.ndimage.gaussian_filter1d(accelero['Norm'].to_numpy(),nfilt)
			Nini=scipy.ndimage.median_filter(accelero['Norm'].to_numpy(),
			filt_acc)
			N=np.copy(Nini)
			# si seuils, modifier N
			print(N.max(),N.min())
			for i in range(len(seuils_acc)-1):
				N[(Nini>seuils_acc[i])&(Nini<seuils_acc[i+1])]=i
		except:
			print('Warning: no accelerometer was read')
			plot_acc=False

	# Read cardio :
	if (plot_bpm)|(plot_alt>0):
		try:
			cardio=read_cardio(patient)
			# filter
			nfilt=50
			#accelero['Norm']=scipy.ndimage.median_filter(accelero['Norm'].to_numpy(),nfilt)
			#N=scipy.ndimage.gaussian_filter1d(accelero['Norm'].to_numpy(),nfilt)
			BPMini=scipy.ndimage.median_filter(cardio['HR (bpm)'].to_numpy(),filt_bpm)
			BPM=np.copy(BPMini)
			#print(np.nanmin(BPM),np.nanmax(BPM))
			for i in range(len(seuils_bpm)-1):
				#print(seuils_bpm[i])
				BPM[(BPMini>seuils_bpm[i])&(BPMini<=seuils_bpm[i+1])]=i
			#print(np.nanmin(BPM),np.nanmax(BPM))
			ALT=scipy.ndimage.median_filter(cardio['Altitude (m)'].to_numpy(),nfilt)*plot_alt
		except:
			print('Warning: no cardio was read')
			plot_bpm=False
			plot_alt=0

	# Intervals
	if not start:
		start = numdate(date_col_num(GluTime_all[:1,:]))
	else:
		start = pd.Timestamp(start)
	if not end:
		end = numdate(date_col_num(GluTime_all[-1:,:]))
	else:
		end = pd.Timestamp(end)
	
	start_1 = pd.Timestamp(start)
	end_1= pd.Timestamp(end)
	
	print(start,end)
	ndays = max(1, int(np.ceil((end_1 - start_1) / pd.Timedelta(days=1))))
	fig,ax=plt.subplots(ndays,1,figsize=(15,ndays*0.2*15),sharex=True)
	ax = np.atleast_1d(ax)
	
	for n in range(ndays):
		day_start = start_1.normalize() + pd.Timedelta(days=n)
		t0 = max(start_1, day_start)
		t1 = min(end_1, day_start + pd.Timedelta(days=1))
		if t1 <= t0:
			ax[n].axis('off')
			continue
# =============================================================================
		# Plot Accelero
# =============================================================================
		if plot_acc:
			isin_pd=(accelero['Date']>=t0)&(accelero['Date']<t1)
			if np.any(isin_pd):
				img=N[isin_pd].reshape(1,-1)
				ax[n].imshow(img,cmap=plt.cm.YlOrBr,extent=[t0,t1,0,superhyper*1.2],alpha=0.4,aspect='auto')


		if plot_bpm:
			isin_pd=(cardio['Time']>=t0)&(cardio['Time']<t1)
			ny=200
			img=BPM[isin_pd].reshape(1,-1)
			if len(seuils_bpm)==0:
				img=BPM[isin_pd].reshape(-1,1)
				img=np.tile(img,ny).T
				x,y=np.meshgrid(np.arange(img.shape[1]),np.arange(img.shape[0]))
				std=ny/800*(BPM[isin_pd]-30).reshape(1,-1)
				img=img*np.exp(-(y-ny/2)**2/(2*std**2))
			#print(BPM[isin_pd])
			#print(n,len(isin_pd))
			ax[n].imshow(img,cmap=plt.cm.PuRd,extent=[t0,t1,0,superhyper*1.2],alpha=0.4,aspect='auto')
	
	
		if plot_alt>0:
			isin_pd=(cardio['Time']>=t0)&(cardio['Time']<t1)
			if np.any(isin_pd):
				ax[n].fill_between(cardio['Time'][isin_pd],np.zeros(sum(isin_pd)),ALT[isin_pd],color='Gray',alpha=0.4,edgecolor='w')
				alt2gly=lambda x:x/plot_alt
				gly2alt=lambda x:x*plot_alt
				secax = ax[n].secondary_yaxis('right', functions=(alt2gly, gly2alt),color='Gray')
				secax.set_xlabel('Altitude [m]')
		#ax[n].xaxis.grid(True, which='minor')
# =============================================================================
# 		Plot Glycemia
# =============================================================================
		# hypo/hyper as areas
		isin_pd=(GluTime_all_pd>=t0)&(GluTime_all_pd<t1)
# 		#ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hypo,np.zeros(sum(isin_pd))+hyper,color='g',alpha=0.1)
# 		# Hyper glycemia
# 		glu_temp=np.copy(GluValue_all[isin_pd,1])
# 		#print(np.sum(glu_temp<hyper),glu_temp)
# 		glu_temp[glu_temp<hyper]=np.nan
# 		ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hyper,glu_temp,color='r',alpha=0.5)
# 		# Hypo glycemia
# 		glu_temp=np.copy(GluValue_all[isin_pd,1])
# 		glu_temp[glu_temp>hypo]=np.nan
		T=[t0,t1]
# 		ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hypo,glu_temp,color='b',alpha=0.5)
		ax[n].plot(T,np.zeros(len(T))+superhyper,'--',linewidth=lw,color='crimson')
		ax[n].plot(T,np.zeros(len(T))+superhypo,'--',linewidth=lw,color='crimson')
		ax[n].plot(T,np.zeros(len(T))+hyper,'-.',color='royalblue',linewidth=lw)
		ax[n].plot(T,np.zeros(len(T))+hypo,'-.',color='royalblue',linewidth=lw)
		ax[n].plot(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),GluValue_all[isin_pd,1],'-',linewidth=lw,color='royalblue')
		ax[n].xaxis.set_major_formatter(xformatter)
		ax[n].set_ylim([GluValue_all[:,1].min(),GluValue_all[:,1].max()])
		ax[n].set_ylim([0,superhyper*1.2])
		ax[n].grid(visible=True)
		ax[n].set_ylabel('Glycemia [mg/dL]',color='royalblue')
		ax[n].tick_params(labelcolor='royalblue')
		ax[n].text(t0,superhyper*1.25,str(GluTime_all_pd[isin_pd][0].strftime('%A %d-%m-%Y')))
		ax[n].set_xlim([t0,t1])

	
	ax[n].set_xlabel('Time')
	if savefig:
		save_folder='./Patients/'+patient+'/'
		if not(os.path.isdir(save_folder)):
			os.mkdir(save_folder)
		tag0 = pd.Timestamp(start_1).strftime('%Y%m%d_%H%M')
		tag1 = pd.Timestamp(end_1).strftime('%Y%m%d_%H%M')
		plt.savefig(save_folder+f'Results_{patient}_{tag0}-{tag1}.pdf',bbox_inches='tight')

def read_glu_pandas(patient,verbose=True):
	
	base_dir = os.path.join(".", "Data", patient)
	if verbose:
		print('Reading Glycemia file in '+base_dir)

	filename,sensors=read_libraries(encoding='utf-8')
# 	print(sensors['delimiter'])
# 	print(filename)
	for i,f in enumerate(filename['Filename']):
		fpath = os.path.join(base_dir, f)
		#print(fpath)
		if os.path.exists(fpath):
			if verbose:
				print('Reading from '+ f)
			sensor=filename['Sensor'][i]
			if verbose:
				print('Sensor format detected:',sensor)
			# Find configuration of sensor in sensors library
			sensor_config=sensors[sensors['Sensor']==sensor]
			#print(sensor_config["delimiter"].iloc[0],sensor_config["skiprows"].iloc[0])
			try:
				Data=pd.read_csv(fpath,delimiter=str(sensor_config["delimiter"].iloc[0]),skiprows=int(sensor_config["skiprows"].iloc[0]),encoding=sensor_config["encoding"].iloc[0],engine='python')
				print(list(Data.columns.values))
				GluValue=pd.to_numeric(Data[sensor_config["col_value"].iloc[0]],downcast='float',errors='coerce').to_numpy()
				GluValue=np.tile(GluValue.reshape(-1,1)*conversion_factor(sensor_config["units"].iloc[0]),2)
				if str(sensor_config["col_datetime"].iloc[0])!='nan':
					TimeTemp=pd.to_datetime(Data[sensor_config["col_datetime"].iloc[0]],errors='coerce',dayfirst=sensor_config["dayfirst"].iloc[0]).fillna(value=pd.Timestamp('20000101'))
				else:
					TimeTemp=pd.to_datetime(Data[sensor_config["col_date"].iloc[0]]+' '+Data[sensor_config["col_time"].iloc[0]],errors='coerce',dayfirst=sensor_config["dayfirst"].iloc[0]).fillna(value=pd.Timestamp('20000101'))
				GluTime=np.vstack((TimeTemp.dt.year.to_numpy(),
						  TimeTemp.dt.month.to_numpy(),
						  TimeTemp.dt.day.to_numpy(),
						  TimeTemp.dt.hour.to_numpy(),
						  TimeTemp.dt.minute.to_numpy(),
						  TimeTemp.dt.second.to_numpy())).T
				GluTime, GluValue=sort_times(GluTime, GluValue)
				GluTime, GluValue=clean_data(GluTime, GluValue)
				print(green+'Reading success ! ', GluTime.shape[0] ,  ' values have been read ! \n'+black)
			except: # Try skiprows+1
				Data=pd.read_csv(fpath,delimiter=str(sensor_config["delimiter"].iloc[0]),skiprows=int(sensor_config["skiprows"].iloc[0]+1),encoding=sensor_config["encoding"].iloc[0],engine='python')
				print(list(Data.columns.values))
				GluValue=pd.to_numeric(Data[sensor_config["col_value"].iloc[0]],downcast='float',errors='coerce').to_numpy()
				GluValue=np.tile(GluValue.reshape(-1,1)*conversion_factor(sensor_config["units"].iloc[0]),2)
				if str(sensor_config["col_datetime"].iloc[0])!='nan':
					TimeTemp=pd.to_datetime(Data[sensor_config["col_datetime"].iloc[0]],errors='coerce',dayfirst=sensor_config["dayfirst"].iloc[0]).fillna(value=pd.Timestamp('20000101'))
				else:
					TimeTemp=pd.to_datetime(Data[sensor_config["col_date"].iloc[0]]+' '+Data[sensor_config["col_time"].iloc[0]],errors='coerce',dayfirst=sensor_config["dayfirst"].iloc[0]).fillna(value=pd.Timestamp('20000101'))
				GluTime=np.vstack((TimeTemp.dt.year.to_numpy(),
						  TimeTemp.dt.month.to_numpy(),
						  TimeTemp.dt.day.to_numpy(),
						  TimeTemp.dt.hour.to_numpy(),
						  TimeTemp.dt.minute.to_numpy(),
						  TimeTemp.dt.second.to_numpy())).T
				GluTime, GluValue=sort_times(GluTime, GluValue)
				GluTime, GluValue=clean_data(GluTime, GluValue)
				print(green+'Reading success ! ', GluTime.shape[0] ,  ' values have been read ! \n'+black)

			return GluTime, GluValue
	print(orange + 'Warning ! No data has been read\n' + black)
	return [],[]



def read_Glu(patient,encoding='utf-8'):
	# General function to read Glycemic data
	try:
		GluTime_all,GluValue_all=read_glu_pandas(patient)
		if len(GluTime_all)==0:
			print('Warning: Reading with read_glu_pandas has failed...')
			print('Trying with read_gluoctave...')
			GluTime_all,GluValue_all=read_glu_octave(patient,encoding=encoding)
	except:
		print('Warning: Reading with read_glu_pandas has failed...')
		print('Trying with read_glu_octave...')
		GluTime_all,GluValue_all=read_glu_octave(patient,encoding=encoding)
	return GluTime_all,GluValue_all

def read_glu_octave(patient,encoding='utf-8'):
	mmolTomg = 18.02  # conversion factor mmol/L to mg/L
	
	GluTime = []
	GluValue = []
	base_dir = os.path.join(".", "Data", patient)
	print('Reading Glycemia file in '+base_dir)
	
	def append_entry(year, month, day, hour, minute, second, val1, val2=None):
		GluTime.append([year, month, day, hour, minute, second])
		if val2 is None:
			GluValue.append([val1, val1])
		else:
			GluValue.append([val1, val2])

	# ===== Format STANDARD =====
	fname="Capteur_standard.csv"
	fpath = os.path.join(base_dir, fname)
	if os.path.exists(fpath):
		print('Reading from '+fname)
		with open(fpath, "r",encoding=FILE_ENCODING['capteur_standard']) as f:
			for line in f:
				line = line.replace(",", ".").strip()
				#print(line)
				try:
					date_str, time_str, val_str = line.split(";")
					dt = datetime.strptime(f"{date_str} {time_str}", "%d/%m/%Y %H:%M:%S")
					append_entry(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, float(val_str))
				except ValueError:
					continue
		GluTime, GluValue=sort_times(GluTime, GluValue)
		return GluTime, GluValue

	# ===== IPRO2 (Capteur.csv + Capteur_ipro2.csv) =====
	for fname in ["Capteur.csv", "Capteur_ipro2.csv"]:
		fpath = os.path.join(base_dir, fname)
		if os.path.exists(fpath):
			with open(fpath, "r",encoding=FILE_ENCODING['capteur_ipro']) as f:
				for line in f:
					
#210	13/05/2014	09:40:19	13/05/2014 09:40	capteur				12,63	62						SGReceivedIPro	"DATE=Tue May 13 09:40:19 UTC 2014;SOURCE=GSR;AMOUNT=62.0;ISIG=13.715428225124402;CALIB_ISIG=12.633458910676769;COMMENT_CODE=NONE;"
					line = line.replace(",", ".")#[1::2]
					# Pattern 1
					m = re.match(
						r"^\d+\s(\d)(\d)/(\d)(\d)/(\d{4})\s(\d)(\d):(\d)(\d):(\d{2})\s\d{2}/\d{2}/\d{4}\s\d{2}:\d{2}\s+capteur\s+([\d.]+)\s+(\d+)",
						line
					)
					if m:
						g = list(map(int, m.groups()[:-2])) + list(map(float, m.groups()[-2:]))
						append_entry(g[4], g[2]*10+g[3], g[0]*10+g[1], g[5]*10+g[6], g[7]*10+g[8], g[9], g[10], g[11])
						continue

					# Pattern 2 (JO1)
					m = re.match(
						r"^\d+\t(\d+)-(\d+)-(\d+)\t(\d+):(\d+):(\d+)\t(\d+)-(\d+)-(\d+)\s+(\d+):(\d+):(\d+)\s+capteur\s+([\d.]+)\s+(\d+)",
						line
					)
					if m:
						g = list(map(float, m.groups()))
						append_entry(int(g[2]) + 2000, int(g[1]), int(g[0]), int(g[3]), int(g[4]), int(g[5]), g[12], g[13])
						continue

					# Pattern 3
					m = re.match(
						r"^\d+\t(\d)(\d)-(\d)(\d)-(\d)(\d)\t(\d)(\d):(\d)(\d):(\d)(\d)\t(\d)(\d)-(\d)(\d)-(\d)(\d)\t(\d)(\d):(\d)(\d):(\d)(\d)\s+capteur\s+([\d.]+)\s+(\d+)",
						line
					)
					if m:
						g = list(map(float, m.groups()))
						append_entry(int(g[4]*10+g[5]) + 2000,
									 int(g[2]*10+g[3]),
									 int(g[0]*10+g[1]),
									 int(g[6]*10+g[7]),
									 int(g[8]*10+g[9]),
									 int(g[10]*10+g[11]),
									 g[24], g[25])
						continue
				GluTime, GluValue=sort_times(GluTime, GluValue)
				return GluTime, GluValue

	# ===== DEXCOM =====
	for fname in ["Capteur.txt", "Capteur_dexcom.txt"]:
		fpath = os.path.join(base_dir, fname)
		if os.path.exists(fpath):
			with open(fpath, "r", encoding=FILE_ENCODING['capteur_dexcom_canada_mmol']) as f:
				for line in f:
					m = re.match(
						r".*?(\d{2})/(\d{2})/(\d{4})\s(\d{2}):(\d{2})\t(\d{2})/(\d{2})/(\d{4})\s(\d{2}):(\d{2})\t(\d+)",
						line
					)
					if m:
						g = list(map(int, m.groups()))
						append_entry(g[7], g[6], g[5], g[8], g[9], 0, g[10])
				GluTime, GluValue=sort_times(GluTime, GluValue)
				return GluTime, GluValue

	# ===== FREESTYLE =====
	fname="Capteur_freestyle.txt"
	fpath = os.path.join(base_dir, fname)
	if os.path.exists(fpath):
		print("Reading from ", fname)
		with open(fpath, "r", encoding=FILE_ENCODING['capteur_freestyle']) as f:
			for line in f:
				# Try several patterns
				patterns = [
					r"^\d+\t(\d{4})/(\d{2})/(\d{2})\s(\d{2}):(\d{2})\t+0\t+(\d+)",
					r"^\d+\s(\d{4})/(\d{2})/(\d{2})\s(\d{2}):(\d{2})\s+0\s+(\d+)",
					r"^\d+\s(\d{2})/(\d{2})/(\d{4})\s(\d{2}):(\d{2})\s+0\s+(\d+)",
					r".{20}(\d{4})/(\d{2})/(\d{2})\s(\d{2}):(\d{2})\s\d+\s(\d+)",
					r"^\d+\s+(\d{2})/(\d{2})/(\d{4})\s(\d{2}):(\d{2})\s+(\d{2})/(\d{2})/(\d{4})\s(\d{2}):(\d{2})\s+0\s+(\d+)"
				]
				matched = False
				for pat in patterns:
					m = re.match(pat, line)
					#print(line,m)
					if m:
						g = list(map(int, m.groups()))
						#print(g)
						if len(g) == 6:
							append_entry(g[0], g[1], g[2], g[3], g[4], 0, g[5])
						elif len(g) == 7:
							append_entry(g[2], g[1], g[0], g[3], g[4], 0, g[6])
						elif len(g) == 5 + 1:
							append_entry(g[0], g[1], g[2], g[3], g[4], 0, g[5])
						elif len(g) == 11:
							append_entry(g[2], g[1], g[0], g[3], g[4], 0, g[10])
						matched = True
						break
			GluTime, GluValue=sort_times(GluTime, GluValue)
			if len(GluTime)==0:
				col_gly='Historique de la glycémie mg/dL'
				col_time="Horodatage de l'appareil"
				skiprows=2
				dayfirst=True
				Data=pd.read_csv(fpath,delimiter='\t',skiprows=skiprows,low_memory=False)#,parse_dates=[['Date', 'Time']])
				print(Data)
				#print(Data)
				GluValue=pd.to_numeric(Data[col_gly],downcast='float',errors='coerce').to_numpy()
				GluValue=np.tile(GluValue.reshape(-1,1),2)
				TimeTemp=pd.to_datetime(Data[col_time],errors='coerce',dayfirst=dayfirst).fillna(value=pd.Timestamp('20000101'))
				GluTime=np.vstack((TimeTemp.dt.year.to_numpy(),
						  TimeTemp.dt.month.to_numpy(),
						  TimeTemp.dt.day.to_numpy(),
						  TimeTemp.dt.hour.to_numpy(),
						  TimeTemp.dt.minute.to_numpy(),
						  TimeTemp.dt.second.to_numpy())).T
			GluTime, GluValue=sort_times(GluTime, GluValue)
			GluTime, GluValue=clean_data(GluTime, GluValue)
			
			print('Successfully Read ',GluTime.shape[0],' values !')
			
			return GluTime, GluValue

	# ===== FREESTYLE PRO =====
	fpath = os.path.join(base_dir, "Capteur_freestyle_pro.csv")
	if os.path.exists(fpath):
		with open(fpath, "r", encoding=FILE_ENCODING['capteur_freestyle_pro']) as f:
			for line in f:
				if len(line) > 35:
					m = re.match(r".{30}(\d{2})-(\d{2})-(\d{4})\s(\d{2}):(\d{2}),\d+,\d+,(\d+)", line)
					if m:
						g = list(map(int, m.groups()))
						append_entry(g[2], g[1], g[0], g[3], g[4], 0, g[5])
						continue
					m = re.match(r"(\d{2})-(\d{2})-(\d{4})\s(\d{2}):(\d{2}),(\d),(\d+),*", line[30:])
					if m:
						g = list(map(int, m.groups()))
						append_entry(g[2], g[1], g[0], g[3], g[4], 0, g[6])
						continue
					m = re.match(r".{34}(\d{2})/(\d{2})/(\d{4})\s(\d{2}):(\d{2});0;(\d+),*", line)
					if m:
						g = list(map(int, m.groups()))
						append_entry(g[2], g[1], g[0], g[3], g[4], 0, g[5])
						continue
			GluTime, GluValue=sort_times(GluTime, GluValue)
			return GluTime, GluValue

	# ===== Medtronics CANADA mmol =====
	fpath = os.path.join(base_dir, "Capteur_medtronics_canada_mmol.csv")
	if os.path.exists(fpath):
		print('Capteur detecté : Capteur_medtronics_canada_mmol.csv')
		with open(fpath, "r", encoding=FILE_ENCODING['capteur_medtronics_canada_mmol']) as f:
			for line in f:
				m = re.match(r"\d+;(\d{2})/(\d{2})/(\d{4});(\d{2}):(\d{2}):(\d{2}).*?([\d.]+);([\d.]+)", line)
				if m:
					g = list(map(float, m.groups()))
					append_entry(int(g[2]), int(g[1]), int(g[0]), int(g[3]), int(g[4]), int(g[5]),
								 g[7], g[6] * mmolTomg)
		if len(GluTime)==0:
			Data=pd.read_csv(fpath,delimiter=';',skiprows=6,low_memory=False)#,parse_dates=[['Date', 'Time']])
			GluValue=pd.to_numeric(Data['Sensor Glucose (mmol/L)'],downcast='float',errors='coerce').to_numpy() * mmolTomg
			GluValue=np.tile(GluValue.reshape(-1,1),2)
			TimeTemp=pd.to_datetime(Data['Date']+' '+Data['Time'],errors='coerce').fillna(value=pd.Timestamp('20000101'))
			GluTime=np.vstack((TimeTemp.dt.year.to_numpy(),
					  TimeTemp.dt.month.to_numpy(),
					  TimeTemp.dt.day.to_numpy(),
					  TimeTemp.dt.hour.to_numpy(),
					  TimeTemp.dt.minute.to_numpy(),
					  TimeTemp.dt.second.to_numpy())).T
		GluTime, GluValue=sort_times(GluTime, GluValue)
		GluTime, GluValue=clean_data(GluTime, GluValue)
		
		print('Successfully Read ',GluTime.shape[0],' values !')
		return GluTime, GluValue

	# ===== Medtronics CANADA mg =====
	fpath = os.path.join(base_dir, "Capteur_medtronics_canada_mg.csv")
	if os.path.exists(fpath):
		with open(fpath, "r", encoding=FILE_ENCODING['capteur_medtronics_canada_mmol']) as f:
			print('Capteur detecté : Capteur_medtronics_canada_mg')
			for line in f:
				line = line.replace(",", ".")
				m = re.match(r"[\d.]+;(\d{4})/(\d{2})/(\d{2});(\d{2}):(\d{2}):(\d{2}).*?([\d.]+);[\d.]+", line)
				if m:
					g = list(map(float, m.groups()))
					append_entry(int(g[0]), int(g[1]), int(g[2]), int(g[3]), int(g[4]), int(g[5]),
								 g[6])
			GluTime, GluValue=sort_times(GluTime, GluValue)
			return GluTime, GluValue

	# ===== Dexcom CANADA mmol =====
	fpath = os.path.join(base_dir, "Capteur_dexcom_canada_mmol.txt")
	if os.path.exists(fpath):
		with open(fpath, "r", encoding=FILE_ENCODING['capteur_dexcom_canada_mmol']) as f:
			print('Capteur detecté : Capteur_dexcom_canada_mmol')
			for line in f:
				#m = re.match(r".*?(\d{4})-(\d{2})-(\d{2})\s(\d{2}):(\d{2}):(\d{2}).*?([\d.]+)", line)
				m = re.match(r"\t\t(\d{4})-(\d{2})-(\d{2})\s(\d{2}):(\d{2}):(\d{2})\t(\d{4})-(\d{2})-(\d{2})\s(\d{2}):(\d{2}):(\d{2})\t(\d.)", line)
				if m:
					g = list(map(float, m.groups()))
					#print(g[6])
					append_entry(int(g[0]), int(g[1]), int(g[2]), int(g[3]), int(g[4]), int(g[5]),
								 g[-1] * mmolTomg)
			GluTime, GluValue=sort_times(GluTime, GluValue)
			return GluTime, GluValue
		
	# ===== Dexcom CANADA mmol =====
	fname="MEDTRONIC_CANADA_EXPORT.csv"
	col_gly='Sensor Glucose (mmol/L)'
	col_date='Date'
	col_time='Time'
	skiprows=6
	fpath = os.path.join(base_dir, fname)
	if os.path.exists(fpath):
		print('Capteur detecté :'+fname)
		Data=pd.read_csv(fpath,delimiter=';',skiprows=skiprows,low_memory=False)#,parse_dates=[['Date', 'Time']])
		print(Data)
		GluValue=pd.to_numeric(Data[col_gly],downcast='float',errors='coerce').to_numpy() * mmolTomg
		GluValue=np.tile(GluValue.reshape(-1,1),2)
		TimeTemp=pd.to_datetime(Data[col_date]+' '+Data[col_time],errors='coerce').fillna(value=pd.Timestamp('20000101'))
		GluTime=np.vstack((TimeTemp.dt.year.to_numpy(),
				  TimeTemp.dt.month.to_numpy(),
				  TimeTemp.dt.day.to_numpy(),
				  TimeTemp.dt.hour.to_numpy(),
				  TimeTemp.dt.minute.to_numpy(),
				  TimeTemp.dt.second.to_numpy())).T
		GluTime, GluValue=sort_times(GluTime, GluValue)
		GluTime, GluValue=clean_data(GluTime, GluValue)
		print('Successfully Read ',GluTime.shape[0],' values !')
		return GluTime, GluValue
	
	# ===== Dexcom CANADA mmol =====
	fname="CLARITY_CANADA_EXPORT.csv"
	col_gly='Glucose Value (mmol/L)'
	col_time='Timestamp (YYYY-MM-DDThh:mm:ss)'
	skiprows=0
	fpath = os.path.join(base_dir, fname)
	if os.path.exists(fpath):
		print('Capteur detecté :'+fname)
		Data=pd.read_csv(fpath,delimiter=';',skiprows=skiprows,low_memory=False)#,parse_dates=[['Date', 'Time']])
		#print(Data)
		GluValue=pd.to_numeric(Data[col_gly],downcast='float',errors='coerce').to_numpy() * mmolTomg
		GluValue=np.tile(GluValue.reshape(-1,1),2)
		TimeTemp=pd.to_datetime(Data[col_time],errors='coerce').fillna(value=pd.Timestamp('20000101'))
		GluTime=np.vstack((TimeTemp.dt.year.to_numpy(),
				  TimeTemp.dt.month.to_numpy(),
				  TimeTemp.dt.day.to_numpy(),
				  TimeTemp.dt.hour.to_numpy(),
				  TimeTemp.dt.minute.to_numpy(),
				  TimeTemp.dt.second.to_numpy())).T
		GluTime, GluValue=sort_times(GluTime, GluValue)
		GluTime, GluValue=clean_data(GluTime, GluValue)
		print('Successfully Read ',GluTime.shape[0],' values !')
		return GluTime, GluValue
	
	# ===== Dexcom CANADA mmol =====
	fname="DIASEND_CANADA_EXPORT.csv"
	col_gly='mmol/L'
	col_time='Time'
	skiprows=4
	fpath = os.path.join(base_dir, fname)
	if os.path.exists(fpath):
		print('Capteur detecté :'+fname)
		Data=pd.read_csv(fpath,delimiter=';',skiprows=skiprows,low_memory=False)#,parse_dates=[['Date', 'Time']])
		#print(Data)
		GluValue=pd.to_numeric(Data[col_gly],downcast='float',errors='coerce').to_numpy() * mmolTomg
		GluValue=np.tile(GluValue.reshape(-1,1),2)
		TimeTemp=pd.to_datetime(Data[col_time],errors='coerce').fillna(value=pd.Timestamp('20000101'))
		GluTime=np.vstack((TimeTemp.dt.year.to_numpy(),
				  TimeTemp.dt.month.to_numpy(),
				  TimeTemp.dt.day.to_numpy(),
				  TimeTemp.dt.hour.to_numpy(),
				  TimeTemp.dt.minute.to_numpy(),
				  TimeTemp.dt.second.to_numpy())).T
		GluTime, GluValue=sort_times(GluTime, GluValue)
		GluTime, GluValue=clean_data(GluTime, GluValue)
		print('Successfully Read ',GluTime.shape[0],' values !')
		return GluTime, GluValue
		
	# ===== TidePool CANADA mmol =====
	fname="TIDEPOOL_CANADA_EXPORT.csv"
	col_gly='Value'
	col_time='Local Time'
	skiprows=0
	fpath = os.path.join(base_dir, fname)
	if os.path.exists(fpath):
		print('Capteur detecté :'+fname)
		Data=pd.read_csv(fpath,delimiter=';',skiprows=skiprows,low_memory=False)#,parse_dates=[['Date', 'Time']])
		#print(Data)
		GluValue=pd.to_numeric(Data[col_gly],downcast='float',errors='coerce').to_numpy() * mmolTomg
		GluValue=np.tile(GluValue.reshape(-1,1),2)
		TimeTemp=pd.to_datetime(Data[col_time],errors='coerce').fillna(value=pd.Timestamp('20000101'))
		GluTime=np.vstack((TimeTemp.dt.year.to_numpy(),
				  TimeTemp.dt.month.to_numpy(),
				  TimeTemp.dt.day.to_numpy(),
				  TimeTemp.dt.hour.to_numpy(),
				  TimeTemp.dt.minute.to_numpy(),
				  TimeTemp.dt.second.to_numpy())).T
		GluTime, GluValue=sort_times(GluTime, GluValue)
		GluTime, GluValue=clean_data(GluTime, GluValue)
		print('Successfully Read ',GluTime.shape[0],' values !')
		return GluTime, GluValue
		
	print('Warning : no glycemia values have been read! \nTry manually converting to standard capteur format :\n'+
	   'capteur_standard.csv : \n%d/%m/%Y;%H:%M:%S;value')
	return [],[]


def read_cardio(patient):
	import glob
	print("cardio")
	files=glob.glob('./Data/'+patient+'/V800_*.CSV')
	for i,filename in enumerate(files):
		#filename='./Data/'+patient+'/V800_12_2021-04-26_09-43-48.CSV'
		# Read first start time
		first = pd.read_csv(filename, nrows=1, sep=',')
		start_datetime=pd.to_datetime(first['Date'].iloc[0]+' '+first['Start time'].iloc[0],dayfirst=True)
		skiprows=2
		separator=','
		# using attribute "usecols" to set the column that will be used because the last column cause problems
		cardio = pd.read_csv(filename, skiprows=skiprows, sep=',')
		cardio['Time']=pd.to_timedelta(cardio['Time'])+start_datetime
		if i==0:
			cardio_all=cardio.copy()
		else:
			#Replace first value by nan tu cut serie
			cardio.loc[0,'HR (bpm)']=np.nan
			cardio.loc[0,'Altitude (m)']=np.nan
			cardio_all=pd.concat([cardio_all,cardio],ignore_index=True)
	return cardio_all.sort_values(by=['Time'])

def read_accelero(patient):
	print("accelero_gt1m method")
	filename='./Data/'+patient+'/Accelero.csv'
	skiprows=10
	separator=','
	# using attribute "usecols" to set the column that will be used because the last column cause problems
	accelero = pd.read_csv(filename, skiprows=skiprows, sep=',',quotechar='#',usecols=["Date"," Time"," Axis1","Axis2","Axis3"])
	
	#accelero = pd.read_csv(filename, skiprows=skiprows, sep='(?<!"),+', usecols=["Date"," Time"," Axis1","Axis2","Axis3"])
	accelero['Date'] = pd.to_datetime(accelero['Date'] + ' '+accelero[' Time'], dayfirst=True, errors='coerce')#"%d/%m/%Y %H:%M:%S
	#print(accelero[' Axis1'].str.replace(r'""', '', regex=True))
	accelero[' Axis1']=pd.to_numeric(accelero[' Axis1'].str.replace(r'"', '', regex=True),errors='coerce')
	accelero['Axis2']=pd.to_numeric(accelero['Axis2'].str.replace(r'"', '', regex=True),errors='coerce')
	accelero['Axis3']=pd.to_numeric(accelero['Axis3'].str.replace(r'"', '', regex=True),errors='coerce')
	accelero['Norm']=np.sqrt(accelero[' Axis1']**2.+accelero['Axis2']**2.+accelero['Axis3']**2.)
	return accelero

def calc_glu(patient='GZ2',
			 GluCible=np.array([[70,140],[70,180],[54,200],[60,300]]),
			 intervals=np.array([]).reshape(-1,2),
			 WRITE='a',
			 encoding='utf-8',
			 verbose=False,
			 interp=False):
# 	patient='GZ2';
# 	GluCible=np.array([[70,140],[70,180],[54,200],[60,300]]);
# 	intervals=np.array([]).reshape(-1,2)
# 	WRITE='a'
# 	encoding='utf-8'
# 	verbose=False
	GluCible=np.array(GluCible)
	# for MAGE
	nstd=1.
	Tav=30.
	# mg/l to mmol
	conv_factor=1.0
	
	# =============================================================================
	# def  calc_Glu_index(patient,interval=[[]],interval_file=False)
	# =============================================================================
	
	#intervals=np.array([["2014-05-12 08:00","2014-05-14 08:00"],["2014-05-12 08:00","2014-05-18 08:00"],["2014-05-12 08:00","2014-05-16 08:00"]])
	
	
	# sensor = pd.read_csv('./data/test/Capteur medtronics.csv',
	# 										 skiprows=11, sep='\t',
	# 										 encoding="utf-16")
	
	# # Turn the date column into datetime objects, format %d-%m-%y %H:%M or %d/%m/%y %H:%M:%s are possible
	# try: sensor['Horodatage'] = pd.to_datetime(sensor['Horodatage'], format='%d/%m/%Y %H:%M')
	# except: sensor['Horodatage'] = pd.to_datetime(sensor['Horodatage'], format='%d-%m-%y %H:%M:%S')
	
	# # Interpolate missing data
	# glycemia_key = 'gly_value'
	# datetime_key = 'Horodatage'
	# glycemia_column = glycemia_key
	# date_column = datetime_key
	# try: sensor[glycemia_key] = sensor['Valeur de capteur (mg/dl)'].interpolate(limit=5)
	# except: sensor[glycemia_key] = sensor['Capteur de glycémie (mg/dl)'].interpolate(limit=5)
	
	# # Raw data
	# Glu=np.array(interval_sensor[glycemia_column])
	# T_i=np.array([datenum(t) for t in interval_sensor[datetime_key]])
	
	# # Interpolation on a rgular time grid (1min)
	# T_interp=np.arange(datenum(interval_sensor[datetime_key].iloc[0]),
	# 									 datenum(interval_sensor[datetime_key].iloc[-1]),1/(60*24)) # reinterpolate data on 1 min intervals
	# Glu_interp=np.interp(T_interp,T_i,Glu)
	
	# To be replaced
	# dataframe=sensor
	# interval_sensor=dataframe
	
	GluTime_all,GluValue_all=read_Glu(patient,encoding=encoding)
	
	if len(GluTime_all)==0:
		return
	
	# Remove some data for testing gaps
# 	GluTime_all=np.delete(GluTime_all,np.arange(20,60),axis=0)
# 	GluValue_all=np.delete(GluValue_all,np.arange(20,60),axis=0)
	
	# Find typical interval bewteen measures
	Time_all=date_col_num(GluTime_all)
	dt_med=np.median(np.diff(Time_all))
	print('Typical time step between measure is ',int(dt_med*24*60),'min')
	# Convert units
	GluValue_all[:, 1] = GluValue_all[:, 1] * conv_factor  # (index 1 because Python is 0-based)
	
	T_interp_all=np.arange(date_col_num(GluTime_all[0,:].reshape(1,-1)),date_col_num(GluTime_all[-1,:].reshape(1,-1)),1/(60*24)) # reinterpolate data on 1 min intervals
	Glu_interp_all=np.interp(T_interp_all,date_col_num(GluTime_all),GluValue_all[:,1])
	Glu_all=np.copy(GluValue_all)
	
	# Date as a panda format
	GluTime_all_pd=numdate(date_col_num(GluTime_all))
	
	dt_med=np.median(np.diff(GluTime_all_pd))
	#plt.hist(np.diff(GluTime_all_pd))
	#print('Typical time step between measure is ',dt_med)
	# Date as a num format
	GluTime_all_num=date_col_num(GluTime_all)
	
	# Tolerance for the absence of measure
	#tlag_Glu = 30 / (60 * 24)  # 30 minutes in days
	tlag_Glu = 1.5 * dt_med / np.timedelta64(1,'D')  # in days
	gaps=np.where(np.diff(GluTime_all_num)>tlag_Glu)[0]
	GluNaN_all=np.ones(Glu_interp_all.shape)
	for g in gaps:
		id_gap=np.where((T_interp_all>GluTime_all_num[g])&(T_interp_all<GluTime_all_num[g+1]))[0]
		GluNaN_all[id_gap]=np.nan

	#Index_all_pd=np.arange(len(GluTime_pd))
	
	# Intervals
	start = numdate(date_col_num(GluTime_all[:1,:]))
	end = numdate(date_col_num(GluTime_all[-1:,:]))
	# add whole period to intervals
	intervals=np.vstack(([start.strftime('%Y-%m-%d %X'),end.strftime('%Y-%m-%d %X')],intervals))
	
	#print(intervals)
#	file_path_inter = Path('Data')/patient/'interval_file.csv'
	file_path_inter = './Data/'+patient+'/interval_file.csv'
	if os.path.exists(file_path_inter):
		inter_df = pd.read_csv(file_path_inter, skiprows=0, sep=';')
		inter_df['Date_start'] = inter_df['Date_start'].astype('string')
		inter_df['Heure_start'] = inter_df['Heure_start'].astype('string')
		#inter_df['timetemps_start'] = pd.to_datetime(inter_df['Date_start']+ ' ' + inter_df['Heure_start'],format='%d/%m/%Y %H:%M:%S',dayfirst=True,errors='coerce')
		inter_df['Date_end'] = inter_df['Date_end'].astype('string')
		inter_df['Heure_end'] = inter_df['Heure_end'].astype('string')
		#inter_df['timetemps_end'] = pd.to_datetime(inter_df['Date_end'] +' ' + inter_df['Heure_end'],format='%d/%m/%Y %H:%M:%S',dayfirst=True,errors='coerce')
		for i in range(len(inter_df)):
			startint=inter_df['Date_start'][i]+' '+ inter_df['Heure_start'][i]
			endint=inter_df['Date_end'][i]+' '+ inter_df['Heure_end'][i]
			intervals=np.vstack((intervals,[startint,endint]))
		
		
		#print("les heures de départ sont :",inter_df['timetemps_start'])
		#print("les heures de fin sont :",inter_df['timetemps_end'])
	else:
		print("Fichier d'interval non trouve")
		
	
	
	IDX_all=[]
	
	print('='*20)
	for interval in intervals:
		isin=(GluTime_all_pd>pd.to_datetime(interval[0]))&(GluTime_all_pd<pd.to_datetime(interval[1]))
		GluTime_pd=GluTime_all_pd[isin]
		Glu_raw=GluValue_all[isin,1]
		#print(GluTime_pd)
		# Loop to remove automesure
		
		
		Index_pd=np.arange(len(GluTime_pd))
		isin_interp=(T_interp_all>datenum(pd.to_datetime(interval[0])))&(T_interp_all<=datenum(pd.to_datetime(interval[1])))
		GluNaN=GluNaN_all[isin_interp]
		Glu_interp=Glu_interp_all[isin_interp]*GluNaN
		T_interp=T_interp_all[isin_interp]
		
		if interp:
			Glu=Glu_interp
		else:
			Glu=Glu_raw
		
		# Dictionnary of index results
		IDX={}
		# =============================================================================
		IDX['start_time'] = interval[0]
		IDX['end_time'] = interval[1]
		print('Interval : ',interval[0],'-',interval[1])
		duration=pd.to_datetime(interval[1])-pd.to_datetime(interval[0])
		
		IDX['duration_hrs'] = duration
		IDX['glycemia_mean'] = np.nanmean(Glu)
		IDX['glycemia_min'] = np.nanmin(Glu)
		IDX['glycemia_max'] = np.nanmax(Glu)
		IDX['glycemia_std'] = np.nanstd(Glu)
		IDX['coef_var'] = (np.nanstd(Glu)/np.nanmean(Glu))*100
		IDX['initial_value'] = Glu[0]
		
		# % of good values
		# old version, computed with interpolated values.
		IDX['%_good_values'] = 100 * (1- np.sum(np.isnan(GluNaN))/len(GluNaN))
		# New version, computed with total interval duration and number of values in interval that have the right dt
		#Do not work because of automesures....
		#nb_values_expected=duration/dt_med
		#nb_values_observed=np.diff(GluTime_pd)>dt_med*0.9
		#print(nb_values_expected)
		#print(len(Glu_raw))
		#IDX['%_good_values'] =  100 * (len(Glu_raw)/nb_values_expected)
		print('Duration: ', duration,'| Good values (%):',int(IDX['%_good_values']))
		# Other index
		nval=Glu.shape[0]
		IDX['glycemia_std_err']=IDX['glycemia_std']/np.sqrt(2*(nval-1)) # Error on std for large nval (web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf)
		IDX['glycemia_mean_err']=IDX['glycemia_mean']/np.sqrt(nval)
		IDX['JIndex']=0.001*(IDX['glycemia_mean']+IDX['glycemia_std'])**2.
		IDX['JIndex_err']=IDX['JIndex']*2*(
		IDX['glycemia_std_err']/IDX['glycemia_std']+IDX['glycemia_mean_err']/IDX['glycemia_mean'])		
		
	# 	Glu_cib_up,Glu_cib_down,Glu_cib_in=[],[],[]
	# 	for i in range(GluCible.shape[0]):
	# 		Glu_cib_up.append(np.nansum(Glu>GluCible[i,1])/np.sum(~np.isnan(Glu))*100) 
	# 		Glu_cib_down.append(np.nansum(Glu<GluCible[i,0])/np.sum(~np.isnan(Glu))*100)
	# 		Glu_cib_in.append(100-Glu_cib_up[i]-Glu_cib_down[i])
	# 	IDX['Glu_cible_up']=Glu_cib_up
	# 	IDX['Glu_cible_down']=Glu_cib_down
	# 	IDX['Glu_cible_in']=Glu_cib_in
	# 	
		#Glu_cib_up,Glu_cib_down,Glu_cib_in=[],[],[]
		for i in range(GluCible.shape[0]):
			glustr='Glu_{:1.0f}-{:1.0f}_'.format(GluCible[i,0],GluCible[i,1])
			IDX[glustr+'up']=np.nansum(Glu>GluCible[i,1])/np.sum(~np.isnan(Glu))*100
			IDX[glustr+'down']=np.nansum(Glu<GluCible[i,0])/np.sum(~np.isnan(Glu))*100
			IDX[glustr+'in']=100-IDX[glustr+'up']-IDX[glustr+'down']
		
		# =============================================================================
		# # Normal hypoglycemia
		# =============================================================================
		#IDX['hypo_Start'],IDX['hypo_Stop']=calc_hypo(Glu_interp,T_interp,[15,15],[54,70])
		#IDX['hypo_Duration']=(IDX['hypo_Stop']-IDX['hypo_Start'])*60*24; # in minutes
		hstart,hstop=calc_hypo(Glu_interp,T_interp,[15,15],[54,70])
		IDX['hypo_Start']=numdate(hstart).round('min')
		IDX['hypo_Stop']=numdate(hstop).round('min')
		IDX['hypo_Duration']=np.uint16((hstop-hstart)*60*24); # in minutes
		if len(IDX['hypo_Duration'])>0:
			IDX['hypo_Duration_mean']=np.mean(IDX['hypo_Duration']) # in minutes
			IDX['hypo_Duration_sum']=np.sum(IDX['hypo_Duration']) # in minutes
			IDX['hypo_Duration_max']=np.max(IDX['hypo_Duration']) # in minutes
			IDX['hypo_Duration_freq']=len(IDX['hypo_Duration'])/(T_interp[-1]-T_interp[0]) # in hypo/day
		else:
			IDX['hypo_Duration_mean']=0 # in minutes
			IDX['hypo_Duration_sum']=0 # in minutes
			IDX['hypo_Duration_max']=0
			IDX['hypo_Duration_freq']=0
		# =============================================================================
		# # Prolongated hypoglycemia
		# =============================================================================
		#IDX['hypo_Start_prol'],IDX['hypo_Stop_prol']=calc_hypo(Glu_interp,T_interp,[120,120],[54,54])
		#IDX['hypo_Duration_prol']=(IDX['hypo_Stop_prol']-IDX['hypo_Start_prol'])*60*24; # in minutes
		hstart,hstop=calc_hypo(Glu_interp,T_interp,[120,120],[54,54])
		IDX['hypo_prol_Start']=numdate(hstart).round('min')
		IDX['hypo_prol_Stop']=numdate(hstop).round('min')
		IDX['hypo_prol_Duration']=np.uint16((hstop-hstart)*60*24); # in minutes
		if len(IDX['hypo_prol_Duration'])>0:
			IDX['hypo_prol_Duration_mean']=np.mean(IDX['hypo_prol_Duration']) # in minutes
			IDX['hypo_prol_Duration_sum']=np.sum(IDX['hypo_prol_Duration']) # in minutes
			IDX['hypo_prol_Duration_max']=np.max(IDX['hypo_prol_Duration']) # in minutes
			IDX['hypo_prol_Duration_freq']=len(IDX['hypo_prol_Duration'])/(T_interp[-1]-T_interp[0]) # in hypo/day
		else:
			IDX['hypo_prol_Duration_mean']=0
			IDX['hypo_prol_Duration_sum']=0
			IDX['hypo_prol_Duration_max']=0
			IDX['hypo_prol_Duration_freq']=0
			
			
		fBG=1.509*(np.log(Glu/conv_factor)**1.084-5.381); # symetric BG (in mg/dl !!)
		rBG=10*fBG**2.; # risk BG
		rlBG=rBG*(fBG<0); # risk BG left
		rhBG=rBG*(fBG>0); # risk BG left
		IDX['LBGI']=np.nanmean(rlBG)
		IDX['LBGI_err']=np.nanstd(rlBG)/np.sqrt(nval)
		IDX['HBGI']=np.nanmean(rhBG)
		IDX['HBGI_err']=np.nanstd(rhBG)/np.sqrt(nval)
		IDX['BGRI']=IDX['LBGI']+IDX['HBGI']
		IDX['BGRI_err']=IDX['LBGI_err']+IDX['HBGI_err']
		
		# GRADE
		GRADE_val=425*(np.log10(np.log10(Glu/18.85))+0.16)**2;
		IDX['GRADE_score_i']=np.nanmean(GRADE_val)
		IDX['GRADE_score_i_err']=np.nanstd(GRADE_val)/np.sqrt(nval);
		HYPO_cont_i=np.nansum(GRADE_val*(Glu<70*conv_factor))/np.nansum(GRADE_val)*100;
		HYPER_cont_i=np.nansum(GRADE_val*(Glu>140*conv_factor))/np.nansum(GRADE_val)*100;
		IDX['EUGLY_cont_i']=100-HYPO_cont_i-HYPER_cont_i;
		
		# M INDEX Schlichtkrull [16] Mirouze et al. [17]
		IVG=np.array([100,120,90]) # for M value
		M_value_IGV_i,M_value_IGV_i_err=[],[]
		for i in range(len(IVG)):
			M_value_IGV_i.append(np.nanmean(np.abs(10*np.log10(Glu/IVG[i])))**3)
			M_value_IGV_i_err.append(np.nanstd(np.abs(10*np.log10(Glu/IVG[i]))**3)/np.sqrt(nval))
		IDX['M_value_IVG']=np.array(M_value_IGV_i)
		IDX['M_value_IVG_err']=np.array(M_value_IGV_i_err)
		
		# Area below curve normalized by the number of correct values
		IDX['Area54']=np.nansum(((54-Glu)>0)*(54-Glu))/np.sum(np.isfinite(Glu))
		IDX['Area70']=np.nansum(((70-Glu)>0)*(70-Glu))/np.sum(np.isfinite(Glu))
		IDX['Area180']=np.nansum(((Glu-180)>0)*(Glu-180))/np.sum(np.isfinite(Glu))
		IDX['Area250']=np.nansum(((Glu-250)>0)*(Glu-250))/np.sum(np.isfinite(Glu))
		
		# ADRR, computed from start_Hour to start_Hour each day
		# Hour of interval is neglected
		
		start_Day=GluTime_pd[0].replace(hour=8,minute=0,second=0)
		end_Day=GluTime_pd[-1].replace(hour=8,minute=0,second=0)
		
		nDays=(end_Day-start_Day).days; # Number of days in interval
		LRi,HRi=[],[]
		for i in range(nDays):
			d0=pd.Timedelta(days=i)
			d1=pd.Timedelta(days=i+1)
		#	print(datetime.datetime.combine(start_Day+d0,start_Hour))
			mask = (GluTime_pd > start_Day+d0) & (GluTime_pd <= start_Day+d1)
			idDays = np.array(Index_pd[mask])#-1 # pd start from 1
			if len(idDays)>0:
				LRi.append(np.nanmax(rlBG[idDays]))
				HRi.append(np.nanmax(rhBG[idDays]))
		LRi=np.array(LRi)
		HRi=np.array(HRi)
		IDX['ADRR']=np.nanmean(LRi+HRi)
		IDX['ADRR_err']=np.nanstd(LRi+HRi)/np.sqrt(nDays)
		
		IDX['CONGAn']=CONGAn(Glu_interp)
		
		Amplitude_real,Time_Amplitude,idAll_filt=calc_MAGE(Glu_interp,T_interp,nstd,Tav)
		IDX['MAGE']=np.nanmean(np.abs(Amplitude_real))
		IDX['MAGE_err']=np.nanstd(np.abs(Amplitude_real))/np.sqrt(len(Amplitude_real))
		IDX['MAGEp']=np.nanmean(Amplitude_real[Amplitude_real>0])
		IDX['MAGEp_err']=np.nanstd(Amplitude_real[Amplitude_real>0])/np.sqrt(np.sum(Amplitude_real>0))
		IDX['MAGEm']=np.nanmean(Amplitude_real[Amplitude_real<0])
		IDX['MAGEm_err']=np.nanstd(Amplitude_real[Amplitude_real<0])/np.sqrt(np.sum(Amplitude_real<0))
		
		# =============================================================================
		# Output
		# =============================================================================
		if verbose:
			for key in list(IDX.keys()):
				try :
					print('{:s} \t \t {:1.3f}'.format(key,IDX[key]))
				except:
					print('{:s} \t \t '.format(key),IDX[key])
		IDX_all.append(IDX)
	
	
	print('='*20)
	# =============================================================================
	# Write CSV result file
	if (WRITE=='a')|(WRITE=='w'):
		sep=';'
		save_folder='./Patients/'+patient+'/'
		if not(os.path.isdir(save_folder)):
			os.mkdir(save_folder)
		with open(save_folder+'result_'+patient+'.csv', WRITE) as f:
			f.write("\n Patient :" + patient + " - calc_Glu.py executed on "+str(pd.Timestamp('now'))+ "\n")
			keys=IDX_all[0].keys()
			for k in keys:
				f.write(k)
				for idx in IDX_all:
					try:
						string=str(idx[k])
						f.write(sep+string.replace('\n',''))
					except:
						f.write('error l1374')
				f.write('\n')
		print('Results written in '+save_folder)
		print('='*20)
	# =============================================================================
	if verbose:
		return IDX_all
	else:
		return


print('GlyApp successfully loaded !')