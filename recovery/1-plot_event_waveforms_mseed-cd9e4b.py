# %% markdown
# # Cut event waveforms from day files
# Searches for events greater than given magnitude and cuts data from day files in SAC format for specified components. Resamples to specified rate, removes instrument response, and saves as SAC file.
# 
# ##### JBR - 2/3/18
# ##### JBR - 2/7/18 : Include option to use GCMT parameters in SAC headers
# %% codecell
%load_ext autoreload
%autoreload
# from setup_parameters_plot import *
import matplotlib.pyplot as plt
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy.io.sac import SACTrace
from obspy.geodetics import gps2dist_azimuth, locations2degrees
from obspy import read
import numpy as np
import pandas as pd
import os
from datetime import datetime
import calendar
from pathlib import Path
# %% codecell

# path2mseed = "/data/irma6/gaherty/youngORCA/OBS/Mseed/" # input miniseed path
# path2mseed = "/data/irma6/gaherty/youngORCA/OBS/" # input miniseed path
path2sac_recov = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/SAC_1Hz_rmresp_killgaps/"
# path2sac_18 = "/Volumes/Russell_2TB/CHAMBO/RESEARCH/PROJ_YoungPacificORCA/DATA/EVENTS/IRIS_XX_5.5_detrend_18sta_Zcorr/"
path2sac_18 = "/Volumes/Russell_2TB/CHAMBO/RESEARCH/PROJ_YoungPacificORCA/DATA/EVENTS/IRIS_XX_5.5_detrend_18sta/"
stalist = "/Volumes/Russell_2TB/YoungORCA_recovered_data/stations_good12.txt"

evname = '201805042233' # M6.9
# evname = '201806212113' # M6.2
# evname = '201806160323' # M6.0
# evname = '201806021153' # M5.9
# evname = '201806180709' # M5.8

amp = 0.3
freqmin=1/150
freqmax=1/20
comp = "BHZ"
network = "XX"

# Recovered stations and hanging stations
recov_stas = ["CC01"]
hang_stas = ['CC02', 'CC04', 'CC11', 'WC03']

# List of differential stations (require multiply by 2 sensitivity)
# DIFF_stas = ["EE03", "CC02", "CC05", "CC11", "WW02"]
# NoClock_stas = ["EC02", "EC04", "EE02", "EE03", "WC01", "WW02"]

# %% codecell
# LOAD STATIONS
inventory = pd.read_csv(stalist, delimiter= '\s+', index_col=False, names=['station','stla','stlo','stel'])

# extract date for naming folders
# date = datetime.strptime(str(tbeg),'%Y-%m-%dT%H:%M:%S.%fZ')
# evname = date.strftime('%Y%m%d%H%M')

# %% codecell
# Loop through good stations 
for ista, station in enumerate(inventory.station):
    inpath = path2sac_18+evname+'/'+evname+'.'+network+'.'+station+'.'+comp+'.sac'
    try:
        if ista == 0:
            st = read(inpath, debug_headers=True)
            npts = st[0].stats.npts
            # Define earthquake parameters
            tbeg = st[0].stats.starttime
            tend = st[0].stats.endtime
            evdp = st[0].stats.sac.evdp/1000
            evla = st[0].stats.sac.evla
            evlo = st[0].stats.sac.evlo
            mag = st[0].stats.sac.mag
        else:
            st += read(inpath, debug_headers=True)
    except Exception:
        print('No file... Skip ' + inpath)
        continue

# Loop through recovered stations
for ista, station in enumerate(recov_stas):
    inpath = '/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/Mseed/CC01/Mseed/CC01.CH2.2018.124.00.21.59.msd'
    inpath_sac = path2sac_recov+station+'/'+station+'.'+str(tbeg.year)+'.'+str(tbeg.julday)+'.00.00.00.'+comp+'.sac'
    
    st_cut = read(inpath)
    temp = read(inpath_sac)
    st_cut[0].stats.sac = temp[0].stats.sac
    
    # t1 = UTCDateTime(tstart)
    # t2 = t1 + 24*60*60        
    st_cut.trim(starttime=tbeg, endtime=tend, pad=True, nearest_sample=False, fill_value=0)
    # Taper new waveform
    st_cut.taper(type="cosine", max_percentage=0.05)
    st_cut.detrend(type='demean')
    st_cut.detrend(type='linear')
    st += st_cut

# %% codecell   
# PLOT!
fl_good = 0
fl_recov = 0
fl_hang = 0
for itr,tr in enumerate(st):    
    # Station info
    stel = tr.stats.sac.stel
    stla = tr.stats.sac.stla
    stlo = tr.stats.sac.stlo    
        
    # station = inventory.station[ista]
    dist, baz, az = gps2dist_azimuth(lat1=stla, lon1=stlo, lat2=evla, lon2=evlo)
    gcarc = locations2degrees(lat1=stla, long1=stlo, lat2=evla, long2=evlo) 
    trfilt = tr.copy()
    trfilt.filter("bandpass", freqmin=freqmin, freqmax=freqmax, zerophase=True)  
    trfilt.normalize()
    
    if itr == 0:
        fig1, ax1 = plt.subplots(figsize=(10,9))
    if tr.stats.station in recov_stas:
        if fl_recov == 0:
            ax1.plot(trfilt.times(),trfilt.data*amp+gcarc,'-r',label='recovered')
            fl_recov = 1
        else:
            ax1.plot(trfilt.times(),trfilt.data*amp+gcarc,'-r')
    elif tr.stats.station in hang_stas:
        if fl_hang == 0:
            ax1.plot(trfilt.times(),trfilt.data*amp+gcarc,'-b',label='hanging')
            fl_hang = 1
        else:
            ax1.plot(trfilt.times(),trfilt.data*amp+gcarc,'-b')
    else:
        if fl_good == 0:
            ax1.plot(trfilt.times(),trfilt.data*amp+gcarc,'-k',label='good')
            fl_good = 1
        else:
            ax1.plot(trfilt.times(),trfilt.data*amp+gcarc,'-k')
    # plt.ylim([0.5e-3, 1])
    plt.xlim([0, 6000])
    #     plt.ylim([-180, -70])
    #             plt.ylim([-200, -80])
    plt.xlabel('Time [s]')
    plt.ylabel('Distance')
    plt.title(str(tbeg)+' ('+str(tbeg.julday)+'): '+comp+' M'+str(mag)+' Depth: '+str(evdp)+' km')
ax1.legend()
plt.show()
