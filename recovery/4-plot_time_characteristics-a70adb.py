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
# %matplotlib inline
# %% codecell

# Load DataFrame
station = 'CC01'
path2pkl = "./pkls/"

df = pd.read_pickle(path2pkl+station+'_data.pkl')
df
# %% codecell
# loop through channels
fig = plt.figure(figsize=(10, 12))
colorday = np.arange(0, 5)
norm = mcolors.Normalize(vmin=colorday.min(), vmax=colorday.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.Set3)
cmap.set_array([])
icount = 0
FS = 13
cmap.to_rgba(4)
for ch, data in df.iteritems():
    ax1 = fig.add_subplot(3, 1, 1)
    ax1.plot(data['jdaynum_mseed_ch'], data['nsec_mseed_actual_ch']/60/60, color=cmap.to_rgba(icount), label=ch, linewidth=4-icount/2)
    # ax1.set_xlabel('Jday', fontsize=FS)
    ax1.set_ylabel('Mseed file length (hours)', fontsize=FS)
    ax1.tick_params(labelsize=FS)
    ax1.grid(True)

    ax2 = fig.add_subplot(3, 1, 2)
    ax2.plot(data['jdaynum_mseed_ch'], data['dt_end_gap_sec_ch'], color=cmap.to_rgba(icount), label=ch, linewidth=4-icount/2)
    # ax2.set_xlabel('Jday', fontsize=FS)
    ax2.set_ylabel('Mseed time file offset (sec)', fontsize=FS)
    ax2.tick_params(labelsize=FS)
    ax2.grid(True)

    ax3 = fig.add_subplot(3, 1, 3)
    ax3.plot(data['jdaynum_mseed_ch'], data['dt_totalshift_sec_save_ch']/60/60, color=cmap.to_rgba(icount), label=ch, linewidth=4-icount/2)
    ax3.set_ylabel('Accumulated gaps (hours)', fontsize=FS)
    ax3.tick_params(labelsize=FS)
    ax3.grid(True)
    ax3.set_xlabel('Jday', fontsize=FS)

    icount += 1
ax1.legend()
plt.show()
# x = df.loc['dt_totalshift_sec_save_ch']['BDH']

fig2 = plt.figure(figsize=(5, 5))
for ch, data in df.iteritems():
    ax1 = fig.add_subplot(1, 1, 1)
    t_mseed = data['tstart_mseed_ch'] + data['npts_mseed_ch'] / data['sr_mseed_ch']
    t_sac = data['tstart_sac_ch'] + data['npts_sac_ch'] / data['sr_sac_ch']
    ax1.plot(data['jdaynum_mseed_ch'], ata['jdaynum_mseed_ch'], color=cmap.to_rgba(icount), label=ch, linewidth=4-icount/2)
    # ax1.set_xlabel('Jday', fontsize=FS)
    ax1.set_ylabel('Mseed file length (hours)', fontsize=FS)
    ax1.tick_params(labelsize=FS)
    ax1.grid(True)

    ax2 = fig.add_subplot(3, 1, 2)
    ax2.plot(data['jdaynum_mseed_ch'], data['dt_end_gap_sec_ch'], color=cmap.to_rgba(icount), label=ch, linewidth=4-icount/2)
    # ax2.set_xlabel('Jday', fontsize=FS)
    ax2.set_ylabel('Mseed time file offset (sec)', fontsize=FS)
    ax2.tick_params(labelsize=FS)
    ax2.grid(True)

    ax3 = fig.add_subplot(3, 1, 3)
    ax3.plot(data['jdaynum_mseed_ch'], data['dt_totalshift_sec_save_ch']/60/60, color=cmap.to_rgba(icount), label=ch, linewidth=4-icount/2)
    ax3.set_ylabel('Accumulated gaps (hours)', fontsize=FS)
    ax3.tick_params(labelsize=FS)
    ax3.grid(True)
    ax3.set_xlabel('Jday', fontsize=FS)

    icount += 1
ax1.legend()
plt.show()