# %% markdown
#### Read single day Mseed files
# Plot a single day of data for each of the four channels

# %% codecell
from obspy import read
from obspy.io.sac import SACTrace
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.colors as mcolors
import os
import glob
import pandas as pd

# %% codecell
# New Mseed files form Sean after filter corrections and clock corrections
path2mseed = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/Mseed/"  # input miniseed path
path2sac = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/SAC_1Hz_rmresp_killgaps/"  # output sac path
path2sta = "/Volumes/Russell_2TB/YoungORCA_recovered_data/stations_recover.txt"  # station file
# path2figdir = "./figs_killgaps/"

# datestr = "2018.124.00.21.59"  # M6.9 starting day
# datestr = "2018.172.00.21.59"  # M6.2 starting day
day = 167  # M6.2 starting day
# datestr = "2018.175.00.21.59"  # ending day
station = "CC01"
chs = ['CH0', 'CH1', 'CH2', 'CH3']

for ich, ch in enumerate(chs):
    for iday, DAY in enumerate([day-1, day, day+1]):
        datestr = "2018."+str(DAY)+".00.21.59"
        mseedfile = path2mseed+station+"/Mseed/"+station+"."+ch+"."+datestr+".msd"
        if iday == 0:
            ch_day = read(mseedfile)
        else:
            ch_day += read(mseedfile)
        # Merge data (interpolate gaps)
        ch_day.merge(method=1, fill_value='interpolate')
    if ich == 0:
        st = ch_day
    else:
        st += ch_day


# %%
st[0].stats.npts
st[0]
# %%
st[1].stats.npts
st[1]
# %%
st[2].stats.npts
st[2]
# %%
st[3].stats.npts
st[3]
# %%
st[0].stats.endtime-st[1].stats.endtime
(st[0].stats.npts-st[1].stats.npts)/st[0].stats.sampling_rate

st.normalize()
st.plot(color='blue')
# %%
st_cut = st.copy()
t1 = st_cut[0].stats.starttime+10*60*60
t2 = t1+3*60*60
st_cut.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)
st_cut.detrend(type='demean')
st_cut.detrend(type='linear')
st_cut.taper(type="cosine", max_percentage=0.05)
st_cut.filter("bandpass", freqmin=1/150, freqmax=1/2, zerophase=True)
st_cut.normalize()
st_cut.plot(color='blue')

