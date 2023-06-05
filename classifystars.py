
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
#from models.barmock import *

AllDiscSNR = read_mock_file(datafile)

classify = False

offset = True

if offset:
    radii = [[5,15],[15,25],[25,35],[35,45]] # the offset bins
    binprefac = 0.1
else:
    radii = [[0,1],[1,2],[2,3],[3,4]]#,[4,5]] # the main bins
    binprefac = 1.0



# print a table of the data
pltkeys = ['f','Lz','Lx', 'Ly','alpha','sxinv', 'syinv', 'szinv']

print('{0:4s}{1:3s}{2:9s}{3:9s} {4:9s} {5:9s} {6:9s} {7:9s} {8:9s} {9:9s}'.format(' ',' ','f','Lz','Lx','Ly','alpha','sx','sy','sz'))

# markdown version for GitHub viewing
if offset:
    f = open(inputdir+'fits/README.md','a')
else:
    f = open(inputdir+'fits/README.md','w')
print('|comp|radii| f | L<sub>z</sub> | L<sub>x</sub> | L<sub>y</sub> | angle | w<sub>z</sub> | w<sub>x</sub> | w<sub>y</sub> |',file=f)
print('|---|---|---| ---| --- | ---| --- | --- | --- | --- |',file=f)

if offset:
    g = open(inputdir+'fits/table.tex','a')
else:
    g = open(inputdir+'fits/table.tex','w')
print('comp &radii & f & $L_z$ & $L_z$& $L_z$ & $\alpha$ & $\sigma_z$ & $\sigma_z$ & $\sigma_z$ \\\\',file=g)

#for irad,rads in enumerate():
for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}".format(modeltag,minrad,maxrad)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)
    print('{0:4d}-{1:1d}'.format(minrad,maxrad))

    for cnum in [0,1,2]:

        # print to markdown
        print('|'+comptag[minrad][cnum],end='|',file=f)
        print('{0:4d}-{1:1d}'.format(minrad,maxrad),end='|',file=f)

        # print to latex
        print(comptag[minrad][cnum],end='&',file=g)
        print('{0:4d}-{1:1d}'.format(minrad,maxrad),end='&',file=g)

        for key in pltkeys:
            if 'inv' in key:
                median = np.sqrt(1./np.nanmedian(CStats[cnum][key]))
                lo = np.sqrt(1./np.nanpercentile(CStats[cnum][key],14.)) - median
                hi = np.sqrt(1./np.nanpercentile(CStats[cnum][key],86.)) - median
                print(' {0:9.4f}'.format(median),end='')
                print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='|',file=f)
                print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='&',file=g)
            elif key=='alpha':
                median = (180./np.pi)*np.nanmedian(CStats[cnum][key])
                lo = (180./np.pi)*np.nanpercentile(CStats[cnum][key],14.) - median
                hi = (180./np.pi)*np.nanpercentile(CStats[cnum][key],86.) - median
                print(' {0:9.4f}'.format(median),end='')
                print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='|',file=f)
                print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='&',file=g)
            elif key=='f':
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median
                print(' {0:9.4f}'.format(median),end='')
                print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,4),np.round(np.abs(lo),4),np.round(np.abs(hi),4)),end='|',file=f)
                print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,4),np.round(np.abs(lo),4),np.round(np.abs(hi),4)),end='&',file=g)
            else: # Lx, Ly, Lz
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median
                print(' {0:9.4f}'.format(median),end='')
                print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='|',file=f)
                print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='&',file=g)
        print(' '+comptag[minrad][cnum])
        print('',file=f)
        print('\\\\',file=g)

f.close()

# generate classifications
if classify:
    # loop through all radii
    for irad,rads in enumerate(radii):

        # get the min/max cylindrical radii in kpc
        minrad,maxrad = rads[0],rads[1]


        # open the correct chain
        directory = inputdir+"fits/{0}_d{1:02d}{2:02d}".format(modeltag,minrad,maxrad)
        inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
        COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

        # define the stars we want to classify
        #criteria = np.where((AllDiscSNR['R']>(minrad)) & (AllDiscSNR['R']<(maxrad)))[0]
        criteria = np.where((AllDiscSNR['R']>(binprefac*minrad)) & (AllDiscSNR['R']<(binprefac*maxrad)))[0]

        # do the probabilistic classification
        # allprobs is the full set of classifications
        # percentileprob is the 50h percentile
        # errorprob is the 1-sigma error on the classification
        allprobs,percentileprob,errorprob = make_all_probabilities(AllDiscSNR,criteria,CStats,nchains=100)

        # print the classification
        disccomp,barcomp,knotcomp = compnum[minrad]
        radii = [minrad,maxrad]
        printdir = inputdir+'classifications/'
        print_classification(AllDiscSNR,criteria,radii,percentileprob,errorprob,disccomp,barcomp,knotcomp,printdir,mockanalysis)
