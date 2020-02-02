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
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.colors as mcolors
import os
import glob
import pandas as pd
# %matplotlib inline
# %% codecell
# Setup paths

# Original Mseed files from Sean before filter corrections and clock corrections
# path2mseed = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/Mseed_old/"  # input miniseed path
# path2sac = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/SAC_1Hz_rmresp_killgaps_old/"  # output sac path
# path2sta = "/Volumes/Russell_2TB/YoungORCA_recovered_data/stations_recover.txt"  # station file
# path2figdir = "./figs_killgaps_old/"

# New Mseed files form Sean after filter corrections and clock corrections
path2mseed = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/Mseed/"  # input miniseed path
path2sac = "/Volumes/Russell_2TB/YoungORCA_recovered_data/DATA/SAC_1Hz_rmresp_killgaps/"  # output sac path
path2sta = "/Volumes/Russell_2TB/YoungORCA_recovered_data/stations_recover.txt"  # station file
path2figdir = "./figs_killgaps/"

# Filter for plotting
freqmax = 1/20
freqmin = 1/150

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
    dt_totalshift_sec_save_ch = {}
    nsec_mseed_actual_ch = {}
    nsec_mseed_apparent_ch = {}
    jdaynum_mseed_ch = {}
    jdaynum_sac_ch = {}
    dt_strt_gap_sec_ch = {}
    dt_end_gap_sec_ch = {}
    mseed_sac_offset_sec_ch = {}
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
        dt_totalshift_sec = 0
        dt_totalshift_sec_save = np.zeros(len(pathlist)-2)
        nsec_mseed_actual = np.zeros(len(pathlist)-2)
        nsec_mseed_apparent = np.zeros(len(pathlist)-2)
        jdaynum_mseed = np.zeros(len(pathlist)-2)
        jdaynum_sac = np.zeros(len(pathlist)-2)
        dt_strt_gap_sec = np.zeros(len(pathlist)-2)
        dt_end_gap_sec = np.zeros(len(pathlist)-2)
        mseed_sac_offset_sec = np.zeros(len(pathlist)-2)
        
        for ifil, path in enumerate(pathlist):
            if (ifil == 0) or (ifil == len(pathlist)-1):
                continue
            iday = ifil-1

            # Load day before, day of, and day after for merging
            st = read(str(pathlist[ifil-1]))
            st += read(str(pathlist[ifil]))
            st += read(str(pathlist[ifil+1]))

            # Check if day already processed
            # if len(glob.glob(path2sac+sta+'/'+sta+'.'+str(st[1].stats.starttime.year)+'.'+'%03i'%(st[1].stats.starttime.julday)+'.*.'+comp+'.sac')) != 0:
            #     print('Skipping... already processed  ' + str(path))
            #     continue
            print("working on: " + str(path))

            # Check that data segments are neighboring
            sr = st[0].stats.sampling_rate
            dt_strt_gap = st[1].stats.starttime-st[0].stats.endtime
            dt_end_gap = st[2].stats.starttime-st[1].stats.endtime
            dt_strt_gap_sec[iday] = dt_strt_gap
            dt_end_gap_sec[iday] = dt_end_gap
            # print(dt_strt_gap)
            # if abs(dt_strt_gap*sr) > 1 or abs(dt_end_gap*sr) > 1:
                # print('Days are not in sequence! ::: '+str(path))
            print('(Day1 strt) - (Day0 end) = ' + str(st[1].stats.starttime-st[0].stats.endtime) + 'sec')
            print('(Day2 strt) - (Day1 end) = ' + str(st[2].stats.starttime-st[1].stats.endtime) + 'sec')
                # continue

            # save original start times!
            otime = st[1].stats.starttime
            oyr = str(st[1].stats.starttime.year)
            ojday = '%03i'%(st[1].stats.starttime.julday)
            ohr = '%02i'%(st[1].stats.starttime.hour)
            omn = '%02i'%(st[1].stats.starttime.minute)
            osec = '%02i'%(st[1].stats.starttime.second)
            nsec_mseed_actual[iday] = st[1].stats.npts / st[1].stats.sampling_rate
            nsec_mseed_apparent[iday] = st[1].stats.endtime - st[1].stats.starttime
            jdaynum_mseed[iday] = ojday
            print('Mseed apparent length: ' + str(nsec_mseed_apparent[iday]) + 's ; Mseed actual length: ' + str(nsec_mseed_actual[iday]))
            # Shift data segments 1 and 2 over
            #
            #           0            1            2
            # Data: ========= <- ========= <- =========
            #                  |            |
            #              dt_strt_gap   dt_end_gap
            #
            st[0].stats.starttime.timestamp = st[0].stats.starttime.timestamp - dt_totalshift_sec
            st[0].stats.endtime.timestamp = st[0].stats.endtime.timestamp - dt_totalshift_sec
            st[1].stats.starttime.timestamp = st[1].stats.starttime.timestamp - dt_strt_gap + 1/sr - dt_totalshift_sec
            st[1].stats.endtime.timestamp = st[1].stats.endtime.timestamp - dt_strt_gap + 1/sr - dt_totalshift_sec
            st[2].stats.starttime.timestamp = st[2].stats.starttime.timestamp - (dt_strt_gap + dt_end_gap) + 2/sr - dt_totalshift_sec
            st[2].stats.endtime.timestamp = st[2].stats.endtime.timestamp - (dt_strt_gap + dt_end_gap) + 2/sr - dt_totalshift_sec

            # Keep track of total time shifted
            if ifil == 1:
                dt_totalshift_sec = dt_totalshift_sec + (dt_strt_gap + dt_end_gap) - 2/sr
            else:  # Make sure not do double-count leading gaps...
                dt_totalshift_sec = dt_totalshift_sec + (dt_end_gap) - 1/sr
            dt_totalshift_sec_save = dt_totalshift_sec

            # Get start time as beginning of middle day
            if ifil == 1:
                refdate = st[1].stats.starttime
                tstart = datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0)
            else:
                tstart = (UTCDateTime(tstart) + 24*60*60).datetime

            # ensure start time in plot equals 0
            dt = st[0].stats.starttime.timestamp

            fig = plt.figure(figsize=(10, 8))
            colorday = np.arange(1, 5)
            norm = mcolors.Normalize(vmin=colorday.min(), vmax=colorday.max())
            cmap = cm.ScalarMappable(norm=norm, cmap=cm.rainbow)
            cmap.set_array([])
            t_zoom = 50
            for i in np.arange(0, 3):
                t = np.linspace(st[i].stats.starttime.timestamp - dt,
                    st[i].stats.endtime.timestamp - dt,
                    st[i].stats.npts)
                ax1 = fig.add_subplot(4, 1, 1)
                ax1.plot(t, st[i].data, color=cmap.to_rgba(i+1))
                ax1.set(ylim=(3*st[i].data.std()+st[i].data.mean(), -3*st[i].data.std()+st[i].data.mean()))
                ax1.set(xlim=(0, 3*24*60*60))
                ax1.set(title=sta+' '+ch+' '+str(st[0].stats.starttime)+' to '+str(st[2].stats.endtime)+': total time shift = '+str(np.round(dt_totalshift_sec/60/60,3))+'hours')
                ax1.set(xticks=[])
                
                ax2 = fig.add_subplot(4, 1, 2)
                tmin = st[0].stats.endtime.timestamp-dt-t_zoom
                tmax = st[1].stats.starttime.timestamp-dt+t_zoom
                tind = (t >= tmin) * (t <= tmax)
                ax2.plot(t[tind], st[i].data[tind], color=cmap.to_rgba(i+1))
                # ax2.set(ylim=(1*st[i].data.std()+st[i].data.mean(), -1*st[i].data.std()+st[i].data.mean()))
                ax2.set(xlim=(tmin, tmax))
                ax2.set(title='start gap='+str(dt_strt_gap)+'s; end gap='+str(dt_end_gap)+'s')
                ax2.set(xticks=[])

                ax3 = fig.add_subplot(4, 1, 3)
                tmin = st[1].stats.endtime.timestamp-dt-t_zoom
                tmax = st[2].stats.starttime.timestamp-dt+t_zoom
                tind = (t >= tmin) * (t <= tmax)
                ax3.plot(t[tind], st[i].data[tind], color=cmap.to_rgba(i+1))
                # ax3.set(ylim=(1*st[i].data.std()+st[i].data.mean(), -1*st[i].data.std()+st[i].data.mean()))
                ax3.set(xlim=(tmin, tmax))
                ax3.set(title='Mseed length: apparent='+str(nsec_mseed_apparent)+'s; actual='+str(nsec_mseed_actual)+'s')
                ax3.set(xticks=[])
            # Merge data (interpolate gaps)
            st.merge(method=1, fill_value='interpolate')

            # Detrend and taper waveform
            st.detrend(type='demean')
            st.detrend(type='linear')

            # make sure trace is 24 hours long
            t1 = UTCDateTime(tstart)
            t2 = t1 + 24*60*60
            st.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)
            
            jdaynum_sac[iday] = st[0].stats.starttime.julday
            mseed_sac_offset_sec[iday] = otime - st[0].stats.starttime
            print('Mseed sac offset: '+str(mseed_sac_offset_sec[iday])+'s')

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

            # Filter the data for plotting
            stfilt = st.copy()
            stfilt[0].filter("bandpass", freqmin=freqmin, freqmax=freqmax, zerophase=True)

            # Plot final trace
            t = np.linspace(stfilt[0].stats.starttime.timestamp - dt,
                stfilt[0].stats.endtime.timestamp - dt,
                stfilt[0].stats.npts)
            ax4 = fig.add_subplot(4, 1, 4)
            ax4.plot(t, stfilt[0].data, color='black')
            ax4.set(ylim=(3*stfilt[0].data.std()+stfilt[0].data.mean(), -3*stfilt[0].data.std()+stfilt[0].data.mean()))
            ax4.set(xlim=(0, 3*24*60*60))
            ax4.set(title=comp+' JDAY:'+str(jdaynum_sac[iday])+'; Mseed sac offset: '+str(mseed_sac_offset_sec[iday]/60/60)+'hr')
            plt.show()
            figdir = path2figdir+sta+'/'
            if 1:
                if not os.path.exists(figdir):
                    os.makedirs(figdir)
                fig.savefig(figdir+sta+'.'+ch+'_'+comp+'.'+oyr+'.'+ojday+'.'+ohr+'.'+omn+'.'+osec+'.pdf')
            fig.clear()
            plt.close(fig)
        dt_totalshift_sec_save_ch[comp] = dt_totalshift_sec_save
        nsec_mseed_actual_ch[comp] = nsec_mseed_actual
        nsec_mseed_apparent_ch[comp] = nsec_mseed_apparent
        jdaynum_mseed_ch[comp] = jdaynum_mseed
        jdaynum_sac_ch[comp] = jdaynum_sac
        dt_strt_gap_sec_ch[comp] = dt_strt_gap_sec
        dt_end_gap_sec_ch[comp] = dt_end_gap_sec
        mseed_sac_offset_sec_ch[comp] = mseed_sac_offset_sec


