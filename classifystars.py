
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


inputdir   = 'data/apogee/'
AllDiscSNR = read_mock_file(inputdir+"APOGEE_all_fehcut_reduceSNR_tagged.txt")

# what is the name of the mock?
mockn = "APOGEE_all_fehcut_reduceSNR"

# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['bar','disc','knot']
comptag[1] = ['disc','bar','knot']
comptag[2] = ['disc','-','bar']
comptag[3] = ['-','-','disc']
comptag[4] = ['-','-','disc']

pltkeys = ['f','Lz','Lx', 'Ly','alpha','sxinv', 'syinv', 'szinv']

print('{0:4s}{1:3s}{2:9s}{3:9s} {4:9s} {5:9s} {6:9s} {7:9s} {8:9s} {9:9s}'.format(' ',' ','f','Lz','Lx','Ly','alpha','sx','sy','sz'))

for minrad in range(0,5):
    maxrad = minrad+1
    directory = "data/apogee/fits/{0}_d{1:02d}{2:02d}".format(mockn,minrad,maxrad)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)
    print('{0:4d}-{1:1d}'.format(minrad,maxrad))
    for cnum in [0,1,2]:
        for key in pltkeys:
            print(' {0:9.4f}'.format(np.nanmedian(CStats[cnum][key])),end='')
        print(' '+comptag[minrad][cnum])
