
import matplotlib.pyplot as plt
# %matplotlib inline
import numpy as np
from astropy.io import fits
import tqdm
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
rc('text')#, usetex=True)

import nmmn.plots
import vorbin

wolfram=nmmn.plots.wolframcmap() # for Mathematica's cmap
parula=nmmn.plots.parulacmap() # for MATLAB's cmap
turbo=nmmn.plots.turbocmap() # Turbo
import warnings
warnings.filterwarnings('ignore')

import pandas as pd

import h5py

from scipy.optimize import curve_fit

from scipy.optimize import fmin, minimize


appendix1="_apog0"
appendix2="_apog1"

# load the csv data
data1 = pd.read_csv('data/barmock/fits/table{}.csv'.format(appendix1))
data2 = pd.read_csv('data/barmock/fits/table{}_consolidated.csv'.format(appendix2))

mask_disc1 = np.where(data1['comp']=='disc')
mask_knot1 = np.where(data1['comp']=='knot')
mask_bar1 = np.where(data1['comp']=='bar')

mask_disc2 = np.where(data2['comp']=='disc')
mask_knot2 = np.where(data2['comp']=='knot')
mask_bar2 = np.where(data2['comp']=='bar')



# this size works well for single-column journal figures
fig = plt.figure(figsize=(3.87,3.0),facecolor='white')

xmin = 0.17
ymin = 0.13
dx = 0.65
dy = 0.83

ax1 = fig.add_axes([xmin+0*dx     ,ymin+0*dy,dx,dy])   # main figure

quantity = 'alpha'
quantity = 'Lz'

ax1.scatter(data1['starmed'].values[mask_bar1],data1[quantity].values[mask_bar1],edgecolor='firebrick',facecolor='None')
ax1.scatter(data2['starmed'].values[mask_bar2],data2[quantity].values[mask_bar2],facecolor='firebrick')

for i in mask_bar1:
    ax1.plot([data1['starmed'].values[i],data1['starmed'].values[i]],[data1[quantity].values[i]+data1[quantity+'-'].values[i],data1[quantity].values[i]+data1[quantity+'+'].values[i]],color='firebrick',linestyle='dashed')

for i in mask_bar2:
    ax1.plot([data2['starmed'].values[i],data2['starmed'].values[i]],[data2[quantity].values[i]+data2[quantity+'-'].values[i],data2[quantity].values[i]+data2[quantity+'+'].values[i]],color='firebrick')

#ax1.axis([0.,1.,0.,1.])
ax1.tick_params(axis="both",direction="in",which="both")
ax1.set_xlabel('radius (kpc)')
ax1.set_ylabel('$\\alpha$')
ax1.set_ylabel('$L_z$')

plt.savefig('/Users/mpetersen/Downloads/mocktest.png',dpi=300)
