


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

mockn = "APOGEE_all_fehcut_reduceSNR"

# what are the (hand-defined) components?
comptag = dict()
comptag[0] = [-1,1,2]#['bar','disc','knot']
comptag[1] = [1,0,2]#['disc','bar','knot']
comptag[2] = [1,-1,0]#['disc','-','bar']
comptag[3] = [-1,-1,1]#['-','-','disc']
comptag[4] = [-1,-1,1]#['-','-','disc']


#         bar    disc   knot
cband = ['blue','black','red']

fig = plt.figure(figsize=(8,3),facecolor='white')
dx = 0.29
dy = 0.70
xmin = 0.1
ymin = 0.15
pointsize = 0.2

ax1 = fig.add_axes([xmin+0*dx,ymin+0*dy,dx,dy])
ax2 = fig.add_axes([xmin+1*dx,ymin+0*dy,dx,dy])
ax3 = fig.add_axes([xmin+2*dx,ymin+0*dy,dx,dy])

ax = [ax1,ax2,ax3]




# this is the angular momentum loop (lz-lx, lz-ly, ly-lx)
for iband,band in enumerate([[3,1],[3,2],[2,1]]):
    print(band[0],band[1])

    # loop through the min/maxrad
    for irad,rads in enumerate([[0,1],[1,2],[2,3],[3,4],[4,5]]):
        minrad,maxrad = rads[0],rads[1]

        directory = "data/apogee/fits/{0}_d{1:02d}{2:02d}".format(mockn,minrad,maxrad)
        inputfile = directory+'/chains/gaussian-post_equal_weights.dat'

        COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

        # loop through the components
        for cnum in [0,1,2]:

            ccnum = comptag[irad][cnum]
            if ccnum < 0: continue

            aellipse,bellipse,color = COMPS[cnum][band[0]+3]**0.5,COMPS[cnum][band[1]+3]**0.5,cband[ccnum]
            xcen,ycen = COMPS[cnum][band[0]],COMPS[cnum][band[1]]

            thellipse1 = np.linspace(0.,np.pi,200)
            xellipse1 = aellipse*np.cos(thellipse1) + xcen
            yellipse1 = bellipse*np.sin(thellipse1) + ycen

            thellipse2 = np.linspace(np.pi,2*np.pi,200)
            xellipse2 = aellipse*np.cos(thellipse2) + xcen
            yellipse2 = bellipse*np.sin(thellipse2) + ycen

            # tilt the ellipse if x/y!
            if iband == 2:
                alpha = COMPS[cnum][7]
                xellipse1tmp = xellipse1*np.cos(alpha) - yellipse1*np.sin(alpha)
                yellipse1tmp = xellipse1*np.sin(alpha) + yellipse1*np.cos(alpha)
                xellipse2tmp = xellipse2*np.cos(alpha) - yellipse2*np.sin(alpha)
                yellipse2tmp = xellipse2*np.sin(alpha) + yellipse2*np.cos(alpha)
                xellipse1 = xellipse1tmp
                yellipse1 = yellipse1tmp
                xellipse2 = xellipse2tmp
                yellipse2 = yellipse2tmp


            #ax[iband].fill_between(xellipse1,yellipse1,yellipse2,edgecolor='none',facecolor=color,alpha=0.5)

            ax[iband].plot(xellipse1,yellipse1,color=color)
            ax[iband].plot(xellipse2,yellipse2,color=color)

            #if iband==2:
            #    ax[iband].text(900,900-(cnum*150),'f={}'.format(np.round(COMPS[cnum][0],3)),color=color,va='top',ha='right')


for ax in [ax1,ax2,ax3]:
    ax.axis([-1000,1000,-1000,1000])

for ax in [ax2,ax3]:
    ax.set_yticklabels(())

ax1.text(-900,900,'Lz-Lx',color='black',va='top',ha='left')
ax2.text(-900,900,'Lz-Ly',color='black',va='top',ha='left')
ax3.text(-900,900,'Ly-Lx',color='black',va='top',ha='left')


#ax2.set_title('{0}<r<{1} kpc ([Fe/H]>-0.4)'.format(minrad,maxrad))



plt.savefig('figures/testfig.png',dpi=300)
