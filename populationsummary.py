

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

# define the solar position relative to the galactic centre
usol,vsol,wsol = 9.,244.,9.
xsol,ysol,zsol = -8.178,0.,0.021

usol,vsol,wsol = 9.,244.,9.
xsol,ysol,zsol = -8.3,0.,0.021

from models.apogee import *
#from models.bulgemock import *
#from models.barmock import *

Stars = read_mock_file(datafile)

classify = True

offset = True

if offset:
    radii = [[5,15],[15,25],[25,35],[35,45]] # the offset bins
    binprefac = 0.1
else:
    radii = [[0,1],[1,2],[2,3],[3,4]]#,[4,5]] # the main bins
    binprefac = 1.0


# what is the mean angular momentum uncertainty in a given bin?

nrealisations = 100

for irad,rads in enumerate(radii):

    # get the min/max cylindrical radii in kpc
    minrad,maxrad = rads[0],rads[1]

    # define the stars we want to classify
    criteria = np.where((Stars['R']>(binprefac*minrad)) & (Stars['R']<(binprefac*maxrad)))[0]

    bigLx = np.zeros([nrealisations,criteria.size])
    bigLy = np.zeros([nrealisations,criteria.size])
    bigLz = np.zeros([nrealisations,criteria.size])

    for n in range(0,nrealisations):
        # make cartesian velocities
        theta  = (np.pi/2.)-Stars['b'][criteria]
        phi    = Stars['l'][criteria]
        vr     = Stars['vlos'][criteria] + np.random.normal(0.,Stars['evlos'][criteria],size=criteria.size)
        newdist = (Stars['dist'][criteria] + np.random.normal(0.,Stars['edist'][criteria],size=criteria.size))
        vphi   = (Stars['pml'][criteria] + np.random.normal(0.,Stars['epml'][criteria],size=criteria.size))*newdist*4.74
        vtheta = -(Stars['pmb'][criteria] + np.random.normal(0.,Stars['epmb'][criteria],size=criteria.size))*newdist*4.74
        v0     = spherical_to_cartesian_velocities(phi,theta,vr,vphi,vtheta)
        x0     = spherical_to_cartesian_positions(phi,theta,newdist)

        x,y,z = x0
        u,v,w = v0
        u += usol
        v += vsol
        w += wsol

        #plt.scatter(Stars['u'],Stars['v'],color='black',s=0.1)

        bigLx[n] = y*w - z*v
        bigLy[n] = z*u - x*w
        bigLz[n] = x*v - y*u

    print('-------------')
    print(minrad)
    print('Lx=',np.nanmedian(np.nanstd(bigLx,axis=0)))
    print('Ly=',np.nanmedian(np.nanstd(bigLy,axis=0)))
    print('Lz=',np.nanmedian(np.nanstd(bigLz,axis=0)))
    #Stars['L'] = np.sqrt(Stars['Lx']*Stars['Lx']+Stars['Ly']*Stars['Ly']+Stars['Lz']*Stars['Lz'])
