# %% markdown
# # Cut event waveforms from day files
# Searches for events greater than given magnitude and cuts data from day files in SAC format for specified components. Resamples to specified rate, removes instrument response, and saves as SAC file.
# 
# ##### JBR - 2/3/18
# ##### JBR - 2/7/18 : Include option to use GCMT parameters in SAC headers
# %% codecell
%load_ext autoreload
%autoreload
from setup_parameters_plot import *
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

# path2mseed = "/data/irma6/gaherty/youngORCA/OBS/Mseed/" # input miniseed path
# path2mseed = "/data/irma6/gaherty/youngORCA/OBS/" # input miniseed path
path2sac_recov = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/SAC_1Hz_rmresp_killgaps/"
path2sac_18 = "/Volumes/Russell_2TB/CHAMBO/RESEARCH/PROJ_YoungPacificORCA/DATA/EVENTS/IRIS_XX_5.5_detrend_18sta/"
stalist = "/Volumes/Russell_2TB/YoungORCA_recovered_data/stations_good12.txt"

evname = "201805042233"
recov_stas = ["CC01"]
network = "XX"

# List of differential stations (require multiply by 2 sensitivity)
DIFF_stas = ["EE03", "CC02", "CC05", "CC11", "WW02"]
# NoClock_stas = ["EC02", "EC04", "EE02", "EE03", "WC01", "WW02"]
# %% codecell
# if not os.path.exists(search_dir):
#     os.makedirs(search_dir)
    
# %% codecell
# LOAD STATIONS
# inventory = client.get_stations(network=network, station=stations, channel=comps[0], starttime=t1, endtime=t2)
# inventory.plot(projection="local",label=False)
# fig = inventory.plot(method="basemap", show=False) 
inventory = pd.read_csv(stalist, delimiter= '\s+', index_col=False, names=['station','stla','stlo','stel'])
# fig.savefig(search_dir+"events.pdf", bbox_inches="tight")

# extract date for naming folders
date = datetime.strptime(str(tbeg),'%Y-%m-%dT%H:%M:%S.%fZ')
evname = date.strftime('%Y%m%d%H%M')

# Loop through stations
for ista, station in enumerate(inventory.station) :
    inpath = path2sac_18+evname+'/'+evname+'.'+network+'.'+station+'.'+comp+'.sac'
    stel = inventory.stel[ista]
    stla = inventory.stla[ista]
    stlo = inventory.stlo[ista]
    
    # Define earthquake parameters
    if ista == 0:
        tbeg = cat.origins[ior].time
        tend = tbeg + trlen
        evdp = cat.origins[ior].depth
        evla = cat.origins[ior].latitude
        evlo = cat.origins[ior].longitude 
        mag = cat.magnitudes[0].mag
    
    # station = inventory.station[ista]
    dist, baz, az = gps2dist_azimuth(lat1=stla, lon1=stlo, lat2=evla, lon2=evlo)
    gcarc = locations2degrees(lat1=stla, long1=stlo, lat2=evla, long2=evlo)

    # Loop through components
    for icomp, ch in enumerate(chs):
        # Get component name (flip BH2 and BH1)
        if ch == "CH0":
            comp = "BH2"
        elif ch == "CH1":
            comp = "BH1"
        elif ch == "CH2":
            comp = "BHZ"
        elif ch == "CH3":
            comp = "BDH"
        # get days before and after event
        datestr_bef = str(cat.origins[ior].time.year)+'.'+'%03i'%((cat.origins[ior].time - 60*60*24).julday)
        datestr = str(cat.origins[ior].time.year)+'.'+'%03i'%(cat.origins[ior].time.julday)
        datestr_aft = str(cat.origins[ior].time.year)+'.'+'%03i'%((cat.origins[ior].time + 60*60*24).julday)
        try:
            if station in NoClock_stas:
                path2mseedsta = path2mseed + 'NoClock/' + station + '/Mseed/'
            else:
                path2mseedsta = path2mseed + 'Mseed/' + station                
    #                 st = client.get_waveforms(network=network, station=station, location="*", channel=comp, starttime=tbeg, endtime=tend, attach_response=True)
            pathlist_bef = sorted(Path(path2mseedsta).glob('**/*'+ch+'*'+datestr_bef+'*.msd'))
            pathlist = sorted(Path(path2mseedsta).glob('**/*'+ch+'*'+datestr+'*.msd'))
            pathlist_aft = sorted(Path(path2mseedsta).glob('**/*'+ch+'*'+datestr_aft+'*.msd'))
            # Load day before, day of, and day after for merging
            st = read(str(pathlist_bef[0]))
            st += read(str(pathlist[0]))
            st += read(str(pathlist_aft[0]))
            # Check that data segments are neighboring
            sr = st[0].stats.sampling_rate
            if (st[1].stats.starttime-st[0].stats.endtime)*sr>50 or (st[2].stats.starttime-st[1].stats.endtime)*sr>50:
                print('::: Days are not in sequence!')
                continue
            # Merge data
            st.merge(method=1)
            st.detrend(type='demean')
            st.detrend(type='linear')
            # Trim to desired length
            st.trim(starttime=tbeg, endtime=tend, pad=True, nearest_sample=False, fill_value=0)
            st.detrend(type='demean')
            st.detrend(type='linear')
            # Taper new waveform
            st.taper(type="cosine",max_percentage=0.05)
        except Exception:
            print('Missing data for station: '+station+' '+comp)
            continue
        if len(st) > 1: # Check for data gaps and fill with 0's
            st.merge(method=1, fill_value=0)
        sr = st[0].stats.sampling_rate
        if is_removeresp:
            try:
                if comp != "BDH":
                    sensitivity = 1/6.85058E-10
                    if inventory.station[ista] in DIFF_stas:
                        sensitivity = sensitivity * 2
        #                 print(STA,' differential... x2 sensitivity')
                    # Displacement paz (m) add extra zero that's 0+0j
                    paz = {
                        'zeros': [0.0 + 0.0j, 0.0 + 0.0j, 0.0 - 0.0j, -1.080000E+02 + 0.000000E+00j, -1.610000E+02 + 0.000000E+00j],
                        'poles': [-1.815000E-02 + 1.799000E-02j, -1.815000E-02 - 1.799000E-02j, -1.730000E+02 + 0.000000E+00j,
                                  -1.960000E+02 + 2.310000E+02j, -1.960000E+02 - 2.310000E+02j, -7.320000E+02 + 1.415000E+03j,
                                  -7.320000E+02 - 1.415000E+03j],
                        'gain': 2.31739E+09,
                        'sensitivity': sensitivity}
                else:
                    # Differential pressure paz (Pa)
                    sensitivity = 1/8.73798E-04
                    paz = {
                        'zeros': [ 0.000000E+00 + 0.000000E+00j],
                        'poles': [-1.256800E-02 + 0.000000E+00j],
                        'gain': 1.00041,
                        'sensitivity': sensitivity}
                st.simulate(paz_remove=paz, paz_simulate=None, pre_filt=[0.001, 0.005, sr/3, sr/2])
                st.detrend(type='demean')
                st.detrend(type='linear')
    #                     st.remove_response(output="DISP", zero_mean=True, taper=True, taper_fraction=0.05, pre_filt=[0.001, 0.005, sr/3, sr/2], water_level=60)
            except Exception:
                print('Failed to remove response: '+evname+' '+station+' '+comp)
                continue
        st.trim(starttime=tbeg, endtime=tend, pad=True, nearest_sample=False, fill_value=0) # make sure correct length
        st.taper(type="cosine",max_percentage=0.05)
        st.detrend(type='demean')
        st.detrend(type='linear')
#             sr_old = st[0].stats.sampling_rate
        if is_downsamp==1:
            st.filter('lowpass', freq=0.4*sr_new, zerophase=True) # anti-alias filter
            st.detrend(type='demean')
            st.detrend(type='linear')
            st.decimate(factor=int(sr/sr_new), no_filter=True) # downsample
            st.detrend(type='demean')
            st.detrend(type='linear')
            st.taper(type="cosine",max_percentage=0.05)
#                 st.resample(sampling_rate=sr_new)

        if ista == 0:
            fig1, ax1 = plt.subplots(figsize=(10,10))
        stcopy = st.copy()
        stcopy[0].filter("bandpass", freqmin=1/150, freqmax=1/20, zerophase=True)  
        stcopy[0].normalize()
        if station in recov_stas:
            ax1.plot(stcopy[0].times(),stcopy[0].data*0.3+gcarc,'-r')
        else:
            ax1.plot(stcopy[0].times(),stcopy[0].data*0.3+gcarc,'-k')
        # plt.ylim([0.5e-3, 1])
        plt.xlim([0, 6000])
    #     plt.ylim([-180, -70])
#             plt.ylim([-200, -80])
        plt.xlabel('Time [s]')
        plt.ylabel('Distance')
        plt.title(comp+' M'+str(mag)+' Depth: '+str(evdp/1000)+' km')
plt.show()

            
# %% codecell
iev
# %% codecell
