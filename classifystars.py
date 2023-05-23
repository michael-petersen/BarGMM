
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

from models.apogee import *
#from models.bulgemock import *

AllDiscSNR = read_mock_file(datafile)

classify = False


# print a table of the data
pltkeys = ['f','Lz','Lx', 'Ly','alpha','sxinv', 'syinv', 'szinv']

print('{0:4s}{1:3s}{2:9s}{3:9s} {4:9s} {5:9s} {6:9s} {7:9s} {8:9s} {9:9s}'.format(' ',' ','f','Lz','Lx','Ly','alpha','sx','sy','sz'))

for minrad in range(0,5):
    maxrad = minrad+1
    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}".format(modeltag,minrad,maxrad)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)
    print('{0:4d}-{1:1d}'.format(minrad,maxrad))
    for cnum in [0,1,2]:
        for key in pltkeys:
            print(' {0:9.4f}'.format(np.nanmedian(CStats[cnum][key])),end='')
        print(' '+comptag[minrad][cnum])



# generate classifications
if classify:
    # loop through all radii
    for irad,rads in enumerate([[0,1],[1,2],[2,3],[3,4],[4,5]]):

        # get the min/max cylindrical radii in kpc
        minrad,maxrad = rads[0],rads[1]

        # open the correct chain
        directory = inputdir+"fits/{0}_d{1:02d}{2:02d}".format(modeltag,minrad,maxrad)
        inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
        COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

        # define the stars we want to classify
        criteria = np.where((AllDiscSNR['R']>(minrad)) & (AllDiscSNR['R']<(maxrad)))[0]

        # do the probabilistic classification
        # allprobs is the full set of classifications
        # percentileprob is the 50h percentile
        # errorprob is the 1-sigma error on the classification
        allprobs,percentileprob,errorprob = make_all_probabilities(AllDiscSNR,criteria,CStats,nchains=100)

        # print the classification
        disccomp,barcomp,knotcomp = compnum[minrad]
        radii = [minrad,maxrad]
        printdir = inputdir+'classifications/'
        print_classification(AllDiscSNR,criteria,radii,percentileprob,errorprob,disccomp,barcomp,knotcomp,printdir)
