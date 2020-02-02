# %% markdown
# ## miniseed to SAC and remove response
#
# Reads miniseed data from the Young Pacific ORCA deployment, chops it into day-long files, and converts it to SAC. This version removes instrument response using user-defined poles and zeros.
#
# Processing steps:
# - Demean
# - Detrend
# - Trim to 24 hour segment
# - Cosine taper
# - Remove instrument response
# - Anti-alias: low pass (corner = 0.4*sr_new Hz)
# - Remove daily fluctuations: high pass (corner = 1/60/60 Hz)
# - Decimate to 1 Hz
# - Convert to SAC and add station lat, lon, depth info
#
#
# To get a left handed system that matches IRIS, need:
# - CH0 -> BH2
# - CH1 -> BH1
# - CH2 -> BHZ
# - CH3 -> BDH
#
# %% codecell
from obspy import read
from obspy.io.sac import SACTrace
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import os
import glob
import pandas as pd
# %% codecell
# Setup paths

# # 12 stations
# path2mseed = "/data/irma6/gaherty/youngORCA/OBS/Mseed/" # input miniseed path
# path2sac = "/data/irma6/jrussel/YoungPacificORCA/SAC_1Hz_18sta_rmresp/" # output sac path
# path2sta = "/Users/jrussel/RESEARCH/PROJ_YoungPacificORCA/DATA/stations_good12.txt" # station file
# isNoClock = 0

# 6 stations with bad clocks
path2mseed = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/Mseed/"  # input miniseed path
path2sac = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/SAC_1Hz_rmresp/"  # output sac path
path2sta = "/Volumes/Russell_2TB/YoungORCA_recovered_data/stations_recover.txt"  # station file

sr_new = 1  # new sample rate in Hz
chs = ["CH0", "CH1", "CH2", "CH3"]

DIFF_stas = ["EE03", "CC02", "CC05", "CC11", "WW02"]
# %% codecell
# Load station file
stala, stalo, staz = np.loadtxt(path2sta, unpack=True, usecols=(1, 2, 3))
if stala.shape == ():  # check for scalar values and make arrays
    stala = np.array([stala]); stalo = np.array([stalo]); staz = np.array([staz])
inventory = pd.read_csv(path2sta, delimiter='\s+', index_col=False, names=['station', 'stla', 'stlo', 'stel'])

# %% codecell
# Loop over stations
for ista,sta in enumerate(inventory.station):
    if not os.path.exists(path2sac+sta):
        os.makedirs(path2sac+sta)
    path2mseedsta = path2mseed + sta + '/Mseed/'
    # Loop over channels
    for ch in chs:
        # Get component name (flip BH2 and BH1)
        if ch == "CH0":
            comp = "BH2"
        elif ch == "CH1":
            comp = "BH1"
        elif ch == "CH2":
            comp = "BHZ"
        elif ch == "CH3":
            comp = "BDH"
        pathlist = sorted(Path(path2mseedsta).glob('**/*'+ch+'*.msd'))

        # Loop over day files and extract days
        for ifil, path in enumerate(pathlist):
            if (ifil == 0) or (ifil == len(pathlist)-1):
                continue

            # Load day before, day of, and day after for merging
            st = read(str(pathlist[ifil-1]))
            st += read(str(pathlist[ifil]))
            st += read(str(pathlist[ifil+1]))

            # Check if day already processed
            if len(glob.glob(path2sac+sta+'/'+sta+'.'+str(st[1].stats.starttime.year)+'.'+'%03i'%(st[1].stats.starttime.julday)+'.*.'+comp+'.sac')) != 0:
                print('Skipping... already processed  ' + str(path))
                continue
            print("working on: " + str(path))

            # Check that data segments are neighboring
            sr = st[0].stats.sampling_rate
            if (st[1].stats.starttime-st[0].stats.endtime)*sr > 50 or (st[2].stats.starttime-st[1].stats.endtime)*sr > 50:
                print('Days are not in sequence! ::: '+str(path))
                print('(Day-1) - (Day) = ' + str((st[1].stats.starttime-st[0].stats.endtime)*sr) + 'sec')
                print('(Day) - (Day+1) = ' + str((st[2].stats.starttime-st[1].stats.endtime)*sr) + 'sec')
                # continue
                quit()

            # Get start time as beginning of middle day
            tstart = st[1]
            refdate = st[1].stats.starttime
            tstart = datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0)

            # Merge data (interpolate gaps)
            st.merge(method=1, fill_value='interpolate')

            # Detrend and taper waveform
            st.detrend(type='demean')
            st.detrend(type='linear')

            # make sure trace is 24 hours long
            t1 = UTCDateTime(tstart)
            t2 = t1 + 24*60*60
            st.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)

            # Taper new waveform
            st.taper(type="cosine", max_percentage=0.05)

            # Remove instrument response using NoMelt
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
                    'gain': 2.316E+09,
                    'sensitivity': sensitivity}
            else:
                # Differential pressure paz (Pa)
                sensitivity = 1/8.73798E-04
                paz = {
                    'zeros': [ 0.000000E+00 + 0.000000E+00j],
                    'poles': [-1.256800E-02 + 0.000000E+00j],
                    'gain': 1.00,
                    'sensitivity': sensitivity}
            st.simulate(paz_remove=paz, paz_simulate=None, pre_filt=[0.001, 0.005, sr/3, sr/2])
#                     st.remove_response(output="DISP", zero_mean=True, taper=True, taper_fraction=0.05, pre_filt=[0.001, 0.005, sr/3, sr/2], water_level=60)

            # Downsample
            st.detrend(type='demean')
            st.detrend(type='linear')
            st.filter('lowpass', freq=0.4*sr_new, zerophase=True)  # anti-alias filter
            st.filter('highpass', freq=1/60/60, zerophase=True)  # Remove daily oscillations
            st.decimate(factor=int(sr/sr_new), no_filter=True)  # downsample
            st.taper(type="cosine", max_percentage=0.05)
            st.detrend(type='demean')
            st.detrend(type='linear')

            # convert to SAC and fill out station/event header info
            sac = SACTrace.from_obspy_trace(st[0])
            sac.stel = staz[ista]
            sac.stla = stala[ista]
            sac.stlo = stalo[ista]
            kcmpnm = comp
            sac.kcmpnm = kcmpnm
            yr = str(st[0].stats.starttime.year)
            jday = '%03i'%(st[0].stats.starttime.julday)
            hr = '%02i'%(st[0].stats.starttime.hour)
            mn = '%02i'%(st[0].stats.starttime.minute)
            sec = '%02i'%(st[0].stats.starttime.second)
            sac_out = path2sac + sta + '/' + sta+'.'+yr+'.'+jday+'.'+hr+'.'+mn+'.'+sec+'.'+kcmpnm+'.sac'
            sac.write(sac_out)

# %% codecell
# sta = "CC04"
# path2tr = path2mseed + sta + "/" + "CC04.CH2.2019.056.00.15.57.msd"
# st = read(path2tr)
# path2tr = path2mseed + sta + "/" + "CC04.CH2.2019.057.00.15.57.msd"
# st += read(path2tr)
# path2tr = path2mseed + sta + "/" + "CC04.CH2.2019.058.00.15.57.msd"
# st += read(path2tr)

# tstart = st[1]
# refdate = st[1].stats.starttime
# tstart = datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0)

# # start time in plot equals 0
# dt = st[0].stats.starttime.timestamp

# ax = plt.subplot(5, 1, 1) # dummy for tying axis
# for i in range(3):
#     plt.subplot(5, 1, i + 1, sharex=ax)
#     t = np.linspace(st[i].stats.starttime.timestamp - dt,
#     st[i].stats.endtime.timestamp - dt,
#     st[i].stats.npts)
#     plt.plot(t, st[i].data)
    
# # Merge the data together and show plot in a similar way
# st.merge(method=1)
# st.detrend(type='demean')
# st.detrend(type='linear')
# plt.subplot(5, 1, 4, sharex=ax)
# t = np.linspace(st[0].stats.starttime.timestamp - dt,
# st[0].stats.endtime.timestamp - dt,
# st[0].stats.npts)
# plt.plot(t, st[0].data, 'r')

# # make sure trace is 24 hours long  
# t1 = UTCDateTime(tstart)
# t2 = t1 + 24*60*60
# st.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)
# # Taper new waveform
# st.taper(type="cosine",max_percentage=0.05)
# # Downsample
# sr_new = 1
# sr = st[0].stats.sampling_rate
# st.filter('lowpass', freq=0.4*sr_new, zerophase=True) # anti-alias filter
# st.filter('highpass', freq=1/60/60, zerophase=True)
# st.decimate(factor=int(sr/sr_new), no_filter=True) # downsample
# plt.subplot(5, 1, 5, sharex=ax)
# t = np.linspace(st[0].stats.starttime.timestamp - dt,
# st[0].stats.endtime.timestamp - dt,
# st[0].stats.npts)
# plt.plot(t, st[0].data, 'k')

# plt.axis('tight')

# plt.show()
