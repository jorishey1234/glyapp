#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 15:25:03 2025

@author: joris
"""

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

def date_col_num(d):
	df = pd.DataFrame(d, columns=['year','month','day','hour','minute','second'])
	num=datenum(df)
	return num

def datenum(d):
	dates = pd.to_datetime(d)
	#res=pd.to_numeric((dates - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s'))/3600/24.
	res1=np.float64((dates - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s'))/3600/24.
	#print(d,'res:',res.to_numpy(),'res1:',res1)
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
	for n in range(1,9): # Hours of shift
		ndt=60*n; #Equivalent in time resolution
		Glu_temp=np.hstack((Glu_interp[ndt:],np.zeros(ndt)+np.nan));
		Dif=Glu_temp-Glu_interp;
		CONGAn.append(np.nanstd(Dif))
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
				 start_Hour=7,hypo=90,hyper=200,
				 plot_acc=True,savefig=True):

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
			nfilt=500
			#accelero['Norm']=scipy.ndimage.median_filter(accelero['Norm'].to_numpy(),nfilt)
			#N=scipy.ndimage.gaussian_filter1d(accelero['Norm'].to_numpy(),nfilt)
			N=scipy.ndimage.median_filter(accelero['Norm'].to_numpy(),nfilt)
		except:
			print('Warning: no accelerometer was read')
			plot_acc=False


	# Intervals
	start = numdate(date_col_num(GluTime_all[:1,:]))
	end = numdate(date_col_num(GluTime_all[-1:,:]))
	
	if start.hour<7:
		start=start+pd.Timedelta(days=-1)
	start_1=start.replace(hour=start_Hour,minute=0,second=0)
	if end.hour>7:
		end=end+pd.Timedelta(days=1)
	end_1=end.replace(hour=start_Hour,minute=0,second=0)
	
	ndays=(end_1-start_1).days
	
	
	fig,ax=plt.subplots(ndays,1,figsize=(15,20),sharex=True)
	for n in range(ndays):
# =============================================================================
		# Plot Accelero
# =============================================================================
		if plot_acc:
			isin_pd=(accelero['Date']>start_1+pd.Timedelta(days=n))&(accelero['Date']<start_1+pd.Timedelta(days=n+1))
			img=N[isin_pd].reshape(1,-1)
			#print(n,len(isin_pd))
			ax[n].imshow(img,cmap=plt.cm.Grays,extent=[start_1,start_1+pd.Timedelta(days=1),0,350],alpha=0.4,aspect='auto')
		#ax[n].xaxis.grid(True, which='minor')
# =============================================================================
# 		Plot Glycemia
# =============================================================================
		isin_pd=(GluTime_all_pd>start_1+pd.Timedelta(days=n))&(GluTime_all_pd<start_1+pd.Timedelta(days=n+1))
		#ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hypo,np.zeros(sum(isin_pd))+hyper,color='g',alpha=0.1)
		# Hyper glycemia
		glu_temp=np.copy(GluValue_all[isin_pd,1])
		#print(np.sum(glu_temp<hyper),glu_temp)
		glu_temp[glu_temp<hyper]=np.nan
		ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hyper,glu_temp,color='r',alpha=0.5)
		# Hypo glycemia
		glu_temp=np.copy(GluValue_all[isin_pd,1])
		glu_temp[glu_temp>hypo]=np.nan
		ax[n].fill_between(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),np.zeros(sum(isin_pd))+hypo,glu_temp,color='b',alpha=0.5)
		
		ax[n].plot(GluTime_all_pd[isin_pd]-pd.Timedelta(days=n),GluValue_all[isin_pd,1],'k-')
		ax[n].xaxis.set_major_formatter(xformatter)
		ax[n].set_ylim([GluValue_all[:,1].min(),GluValue_all[:,1].max()])
		ax[n].set_ylim([0,350])
		ax[n].grid(visible=True)
		ax[n].text(start_1,350,str(GluTime_all_pd[isin_pd][0].strftime('%A %d-%m-%Y')))
	ax[0].set_xlim([start_1,start_1+pd.Timedelta(days=1)])
	
	ax[n].set_xlabel('Time')
	if savefig:
		save_folder='./Patients/'+patient+'/'
		if not(os.path.isdir(save_folder)):
			os.mkdir(save_folder)
		plt.savefig(save_folder+'Results_'+patient+'.pdf',bbox_inches='tight')


def read_Glu(patient,encoding='utf-8'):
	
	
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

	# ===== Format STANDARD =====
	fpath = os.path.join(base_dir, "Capteur_standard.csv")
	if os.path.exists(fpath):
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
	fpath = os.path.join(base_dir, "Capteur_freestyle.txt")
	if os.path.exists(fpath):
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
		with open(fpath, "r", encoding=FILE_ENCODING['capteur_medtronics_canada_mmol']) as f:
			for line in f:
				m = re.match(r"\d+;(\d{2})/(\d{2})/(\d{4});(\d{2}):(\d{2}):(\d{2}).*?([\d.]+);([\d.]+)", line)
				if m:
					g = list(map(float, m.groups()))
					append_entry(int(g[2]), int(g[1]), int(g[0]), int(g[3]), int(g[4]), int(g[5]),
								 g[7], g[6] * mmolTomg)
			GluTime, GluValue=sort_times(GluTime, GluValue)
			return GluTime, GluValue

	# ===== Medtronics CANADA mg =====
	fpath = os.path.join(base_dir, "Capteur_medtronics_canada_mg.csv")
	if os.path.exists(fpath):
		with open(fpath, "r", encoding=FILE_ENCODING['capteur_medtronics_canada_mmol']) as f:
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
			for line in f:
				#m = re.match(r".*?(\d{4})-(\d{2})-(\d{2})\s(\d{2}):(\d{2}):(\d{2}).*?([\d.]+)", line)
				m = re.match(r"\t\t(\d{4})-(\d{2})-(\d{2})\s(\d{2}):(\d{2}):(\d{2})\t(\d{4})-(\d{2})-(\d{2})\s(\d{2}):(\d{2}):(\d{2})\t(\d.)", line)
				if m:
					g = list(map(float, m.groups()))
					print(g[6])
					append_entry(int(g[0]), int(g[1]), int(g[2]), int(g[3]), int(g[4]), int(g[5]),
								 g[-1] * mmolTomg)
			GluTime, GluValue=sort_times(GluTime, GluValue)
			return GluTime, GluValue
	print('Warning : no glycemia values have been read! \nTry manually converting to standard capteur format :\n'+
	   'capteur_standard.csv : \n%d/%m/%Y;%H:%M:%S;value')
	return [],[]


def read_accelero(patient):
	print("accelero_gt1m method")
	filename='./Data/'+patient+'/Accelero.csv'
	skiprows=10
	separator=','
	# using attribute "usecols" to set the column that will be used because the last column cause problems
	accelero = pd.read_csv(filename, skiprows=skiprows, sep=',',quotechar='#',usecols=["Date"," Time"," Axis1","Axis2","Axis3"])
	
	#accelero = pd.read_csv(filename, skiprows=skiprows, sep='(?<!"),+', usecols=["Date"," Time"," Axis1","Axis2","Axis3"])
	accelero['Date'] = pd.to_datetime(accelero['Date'] + ' '+accelero[' Time'], format='"%d/%m/%Y %H:%M:%S')
	
	accelero[' Axis1']=pd.to_numeric(accelero[' Axis1'].str.replace(r'""', '', regex=True))
	accelero['Axis2']=pd.to_numeric(accelero['Axis2'].str.replace(r'""', '', regex=True))
	accelero['Axis3']=pd.to_numeric(accelero['Axis3'].str.replace(r'""', '', regex=True))
	accelero['Norm']=np.sqrt(accelero[' Axis1']**2.+accelero['Axis2']**2.+accelero['Axis3']**2.)
	return accelero

def calc_glu(patient='GZ2',
			 GluCible=np.array([[70,140],[70,180],[54,200],[60,300]]),
			 intervals=np.array([]).reshape(-1,2),
			 WRITE='a',
			 encoding='utf-8'):
	
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
	
	T_interp_all=np.arange(date_col_num(GluTime_all[0,:].reshape(1,-1)),date_col_num(GluTime_all[-1,:].reshape(1,-1)),1/(60*24)) # reinterpolate data on 1 min intervals
	Glu_interp_all=np.interp(T_interp_all,date_col_num(GluTime_all),GluValue_all[:,1])
	Glu_all=np.copy(GluValue_all)
	
	# Date as a panda format
	GluTime_all_pd=numdate(date_col_num(GluTime_all))
	#Index_all_pd=np.arange(len(GluTime_pd))
	
	# Intervals
	start = numdate(date_col_num(GluTime_all[:1,:]))
	end = numdate(date_col_num(GluTime_all[-1:,:]))
	# add whole period to intervals
	intervals=np.vstack(([start.strftime('%Y-%m-%d %X'),end.strftime('%Y-%m-%d %X')],intervals))
	
	
	IDX_all=[]
	
	for interval in intervals:
		isin=(GluTime_all_pd>pd.to_datetime(interval[0]))&(GluTime_all_pd<=pd.to_datetime(interval[1]))
		GluTime_pd=GluTime_all_pd[isin]
		Glu=GluValue_all[isin,1]
		Index_pd=np.arange(len(GluTime_pd))
		isin_interp=(T_interp_all>datenum(pd.to_datetime(interval[0])))&(T_interp_all<=datenum(pd.to_datetime(interval[1])))
		Glu_interp=Glu_interp_all[isin_interp]
		T_interp=T_interp_all[isin_interp]
		
		# Dictionnary of index results
		IDX={}
		# =============================================================================
		IDX['start_time'] = interval[0]
		IDX['end_time'] = interval[1]
		IDX['glycemia_mean'] = np.nanmean(Glu)
		IDX['glycemia_min'] = np.nanmin(Glu)
		IDX['glycemia_max'] = np.nanmax(Glu)
		IDX['glycemia_std'] = np.nanstd(Glu)
		
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
		IDX['hypo_Start'],IDX['hypo_Stop']=calc_hypo(Glu_interp,T_interp,[15,15],[54,70])
		IDX['hypo_Duration']=(IDX['hypo_Stop']-IDX['hypo_Start'])*60*24; # in minutes
		
		# =============================================================================
		# # Prolongated hypoglycemia
		# =============================================================================
		IDX['hypo_Start_prol'],IDX['hypo_Stop_prol']=calc_hypo(Glu_interp,T_interp,[120,120],[54,54])
		IDX['hypo_Duration_prol']=(IDX['hypo_Stop_prol']-IDX['hypo_Start_prol'])*60*24; # in minutes
		
		
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
		for key in list(IDX.keys()):
			try :
				print('{:s} \t \t {:1.3f}'.format(key,IDX[key]))
			except:
				print('{:s} \t \t '.format(key),IDX[key])
		IDX_all.append(IDX)
	
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
					string=str(idx[k])
					f.write(sep+string.replace('\n',''))
				f.write('\n')
	# =============================================================================
	return IDX_all
