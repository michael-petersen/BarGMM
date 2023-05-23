
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
#from models.bulgemock import *

AllDiscSNR = read_mock_file(datafile)

# set the plotting threshold
threshold = 0.8

for component in ['bar','disc','knot']:

    fig = plt.figure(figsize=(6,6),facecolor='white')

    xmin = 0.1
    ymin = 0.1
    dx = 0.43
    dy = 0.43
    ax1 = fig.add_axes([xmin+0*dx,ymin+0*dy,dx,dy])
    ax2 = fig.add_axes([xmin+1*dx,ymin+0*dy,dx,dy])
    ax3 = fig.add_axes([xmin+0*dx,ymin+1*dy,dx,dy])

    #
    # make a plot of the positions
    #

    # loop through all radii
    for irad,rads in enumerate([[0,1],[1,2],[2,3],[3,4],[4,5]]):

        # get the min/max cylindrical radii in kpc
        minrad,maxrad = rads[0],rads[1]

        # define the stars we want the classifications for
        criteria = np.where((AllDiscSNR['R']>(minrad)) & (AllDiscSNR['R']<(maxrad)))[0]

        # read the classifications
        A = np.genfromtxt(inputdir+'classifications/3Component_AllFeHCutMembership_Percentiles_reduceSNR_r{}R{}_cyl.csv'.format(minrad,maxrad),skip_header=1,delimiter=',')

        # this tracks how many stars had failed classifications
        # failed classifications means they were assigned to nonsense components that were not tracked
        print('classification size',A[:,0].size)
        print('criteria size',AllDiscSNR['l'][criteria].size)

        # plot all data as background
        ax1.scatter(A[:,7],A[:,8],color='grey',s=0.5,zorder=-99)
        ax2.scatter(A[:,9],A[:,8],color='grey',s=0.5,zorder=-99)
        ax3.scatter(A[:,7],A[:,9],color='grey',s=0.5,zorder=-99)

        if component == 'disc':
            clr = 'black'
            indx = 5
        elif component == 'knot':
            clr = 'red'
            indx = 1
        elif component == 'bar':
            clr = 'blue'
            indx = 3
        else:
            print('unrecognised component')
            break

        Scriteria = A[:,indx]>threshold

        ax1.scatter(A[:,7][Scriteria],A[:,8][Scriteria],color=clr,s=1.)
        ax2.scatter(A[:,9][Scriteria],A[:,8][Scriteria],color=clr,s=1.)
        ax3.scatter(A[:,7][Scriteria],A[:,9][Scriteria],color=clr,s=1.)


    for ax in [ax1,ax2,ax3]:
        ax.axis([-5.,5,-5.,5.])
        ax.tick_params(axis="both",direction="in",top=True,right=True)

    ax2.set_yticklabels(())
    ax3.set_xticklabels(())

    ax1.set_xlabel('X (kpc)')
    ax1.set_ylabel('Y (kpc)')
    ax2.set_xlabel('Z (kpc)')
    ax3.set_ylabel('Z (kpc)')

    plt.savefig('figures/allposview_{}_{}.png'.format(component,modelname),dpi=300)

    #
    # make a plot of the positions
    #

    fig = plt.figure(figsize=(6.5,6),facecolor='white')

    xmin = 0.13
    ymin = 0.1
    dx = 0.43
    dy = 0.43
    ax1 = fig.add_axes([xmin+0*dx,ymin+0*dy,dx,dy]) # lower left
    ax2 = fig.add_axes([xmin+1*dx,ymin+0*dy,dx,dy]) # lower right
    ax3 = fig.add_axes([xmin+0*dx,ymin+1*dy,dx,dy]) # upper left


    # loop through all radii
    for irad,rads in enumerate([[0,1],[1,2],[2,3],[3,4],[4,5]]):

        # get the min/max cylindrical radii in kpc
        minrad,maxrad = rads[0],rads[1]

        # define the stars we want the classifications for
        criteria = np.where((AllDiscSNR['R']>(minrad)) & (AllDiscSNR['R']<(maxrad)))[0]

        # read the classifications
        A = np.genfromtxt(inputdir+'classifications/3Component_AllFeHCutMembership_Percentiles_reduceSNR_r{}R{}_cyl.csv'.format(minrad,maxrad),skip_header=1,delimiter=',')

        # this tracks how many stars had failed classifications
        # failed classifications means they were assigned to nonsense components that were not tracked
        print('classification size',A[:,0].size)
        print('criteria size',AllDiscSNR['l'][criteria].size)

        # plot all data as background
        ax1.scatter(A[:,10],A[:,11],color='grey',s=0.5,zorder=-99) # Lx-Ly
        ax2.scatter(A[:,12],A[:,11],color='grey',s=0.5,zorder=-99) # Lz-Ly
        ax3.scatter(A[:,10],A[:,12],color='grey',s=0.5,zorder=-99) # Lx-Lz

        if component == 'disc':
            clr = 'black'
            indx = 5
        elif component == 'knot':
            clr = 'red'
            indx = 1
        elif component == 'bar':
            clr = 'blue'
            indx = 3
        else:
            print('unrecognised component')
            break

        Scriteria = A[:,indx]>threshold

        ax1.scatter(A[:,10][Scriteria],A[:,11][Scriteria],color=clr,s=1.)
        ax2.scatter(A[:,12][Scriteria],A[:,11][Scriteria],color=clr,s=1.)
        ax3.scatter(A[:,10][Scriteria],A[:,12][Scriteria],color=clr,s=1.)


    for ax in [ax1,ax2,ax3]:
        ax.axis([-1000.,1000.,-1000.,1000.])
        ax.tick_params(axis="both",direction="in",top=True,right=True)

    ax2.set_yticklabels(())
    ax3.set_xticklabels(())

    ax1.set_xlabel('Lx (kpc km/s)')
    ax1.set_ylabel('Ly (kpc km/s)')
    ax2.set_xlabel('Lz (kpc km/s)')
    ax3.set_ylabel('Lz (kpc km/s)')

    plt.savefig('figures/allangmomview_{}_{}.png'.format(component,modelname),dpi=300)
