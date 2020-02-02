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

datestr = "2018.124.00.21.59"  # starting day
# datestr = "2018.175.00.21.59"  # ending day
station = "CC01"

# Miniseed paths
mseedfile_ch0 = path2mseed+station+"/Mseed/"+station+".CH0."+datestr+".msd"
mseedfile_ch1 = path2mseed+station+"/Mseed/"+station+".CH1."+datestr+".msd"
mseedfile_ch2 = path2mseed+station+"/Mseed/"+station+".CH2."+datestr+".msd"
mseedfile_ch3 = path2mseed+station+"/Mseed/"+station+".CH3."+datestr+".msd"

# Read miniseed files
st = read(mseedfile_ch0)
st += read(mseedfile_ch1)
st += read(mseedfile_ch2)
st += read(mseedfile_ch3)

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
st.plot()
# %%
st_cut = st.copy()
t1 = st_cut[0].stats.starttime+10.5*60*60
t2 = t1+2*60*60
st_cut.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)
st_cut.filter("bandpass", freqmin=1/70, freqmax=1/20, zerophase=True)
st_cut.plot()

