
# basic imports
import numpy as np
from numpy.linalg import eig, inv

# plotting elements
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
mpl.rcParams['font.weight'] = 'medium';mpl.rcParams['xtick.labelsize'] = 10;mpl.rcParams['ytick.labelsize'] = 10

# exptool imports
from exptool.utils import kde_3d
from exptool.analysis import pattern
from exptool.observables import transform
from exptool.io import particle

import scipy.interpolate as interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import scipy

from src.localtools import *
from src.fitmanagement import *

# read in the model that you want
from models.apogee import *
#from models.barmock import *

import h5py

Stars = read_mock_file(datafile)

f = h5py.File("data/apogee/classifications/AllClassifications_{0}{1}.h5".format(modeltag,appendix),"r")
g = h5py.File("data/apogee/classifications/AllClassifications_{0}{1}_500.h5".format(modeltag,appendix),"r")

keylst = [i for i in g.keys()]
print(keylst[0:15])

plt.figure(figsize=(4,3))

"""
#for key in f.keys():
for ikey,key in enumerate(keylst[::800]):
    clr = cm.viridis((800*ikey)/len(keylst),1.)
    print(f[key][:,0],f[key][0,3],g[key][0,3])
    if '*' not in key:
        key2 = key+'*'
    else:
        key2 = key.rstrip('*')
    for i in range(20,500,40):
        plt.scatter(i,np.nanmedian(g[key][0:i,0]),facecolor=clr,edgecolor='None',s=3)
        try:
            plt.scatter(i,np.nanmedian(g[key2][0:i,0]),facecolor='None',edgecolor=clr,s=5)
        except:
            print('no counterpart',key2)       #plt.scatter(i,np.nanpercentile(g[key][0:i,0],16.),color='black',s=3)
"""

# compare the two different classification schemes
for ikey,key in enumerate(keylst[::200]):
    clr = cm.viridis((200*ikey)/len(keylst),1.)
    print(f[key][:,0],f[key][0,3],g[key][0,3])
    if '*' not in key:
        key2 = key+'*'
    else:
        key2 = key.rstrip('*')
    try:
        plt.scatter(np.nanmedian(g[key][:,0]),np.nanmedian(g[key2][:,0]),facecolor='None',edgecolor=clr,s=3,lw=1.)
    except:
        pass


plt.savefig('/Users/mpetersen/Downloads/medconvergence.png',dpi=300)
