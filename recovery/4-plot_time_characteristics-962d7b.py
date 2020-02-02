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

# %% codecell