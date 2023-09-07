
# basic imports
import numpy as np
from numpy.linalg import eig, inv

# printing imports
import h5py

# plotting elements
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
mpl.rcParams['font.weight'] = 'medium';mpl.rcParams['xtick.labelsize'] = 10;mpl.rcParams['ytick.labelsize'] = 10
plt.ion()

from src.localtools import *
from src.fitmanagement import *

from models.apogee import *
from models.bulgemock import *
#from models.barmock import *

Stars = read_mock_file(datafile)


classify = True


# print a table of the data
pltkeys = ['f','Lx', 'Ly','Lz','alpha','sxinv', 'syinv', 'szinv']


fig = plt.figure(figsize=(12.5,3.5),facecolor='white')

fig = plt.gcf()
xmin = 0.04
ymin = 0.12
dx = 0.18
dy = 0.42
xbuf = 0.06

# upper row
ax1 = fig.add_axes([xmin+0*(dx+xbuf),ymin+1*dy,dx,dy]) # f
ax2 = fig.add_axes([xmin+1*(dx+xbuf),ymin+1*dy,dx,dy]) # Lx
ax3 = fig.add_axes([xmin+2*(dx+xbuf),ymin+1*dy,dx,dy]) # Ly
ax4 = fig.add_axes([xmin+3*(dx+xbuf),ymin+1*dy,dx,dy]) # Lz

ax5 = fig.add_axes([xmin+0*(dx+xbuf),ymin+0*dy,dx,dy]) # alpha
ax6 = fig.add_axes([xmin+1*(dx+xbuf),ymin+0*dy,dx,dy]) # sigmax
ax7 = fig.add_axes([xmin+2*(dx+xbuf),ymin+0*dy,dx,dy]) # sigmay
ax8 = fig.add_axes([xmin+3*(dx+xbuf),ymin+0*dy,dx,dy]) # sigmaz


axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]


cwheel = ['blue','black','red']


#for irad,rads in enumerate():
for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

    criteria = np.where((Stars['R']>(binprefacs[irad]*minrad)) & (Stars['R']<(binprefacs[irad]*maxrad)))[0]

    for cnum in [0,1,2]:


        for ikey,key in enumerate(pltkeys):
            if 'inv' in key:
                median = np.sqrt(1./np.nanmedian(CStats[cnum][key]))
                lo = np.sqrt(1./np.nanpercentile(CStats[cnum][key],14.)) - median
                hi = np.sqrt(1./np.nanpercentile(CStats[cnum][key],86.)) - median

            elif key=='alpha':
                median = (180./np.pi)*np.nanmedian(CStats[cnum][key])
                lo = (180./np.pi)*np.nanpercentile(CStats[cnum][key],14.) - median
                hi = (180./np.pi)*np.nanpercentile(CStats[cnum][key],86.) - median
            elif key=='f':
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median

            else: # Lx, Ly, Lz
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median

            radval = np.nanmedian(Stars['R'][criteria])
            if comptag[minrad][cnum]=='bar':
                axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='black')
                axlist[ikey].scatter(radval,median,marker='X',edgecolor='none',facecolor='black')
            elif comptag[minrad][cnum]=='disc':
                axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='blue')
                axlist[ikey].scatter(radval,median,marker='X',edgecolor='none',facecolor='blue')
            elif comptag[minrad][cnum]=='knot':
                axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='red')
                axlist[ikey].scatter(radval,median,marker='X',edgecolor='none',facecolor='red')


appendix="_apog0"
# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['disc','bar','knot']
comptag[1] = ['disc','bar','knot']
comptag[2] = ['disc','-','bar']
comptag[3] = ['-','disc','bar']
comptag[4] = ['-','-','disc']

# by component number, which component is [disc,bar,knot]?
compnum = dict()
compnum[0] = [0,1,2]
compnum[1] = [0,1,2]
compnum[2] = [0,2,-1]
compnum[3] = [1,2,-1]
compnum[4] = [2,-1,-1]

# for ellipse tracing
# by [bar,disc,knot], which component is which?
complist = dict()
complist[0] = [1,0,2]
complist[1] = [1,0,2]
complist[2] = [1,-1,0]
complist[3] = [-1,1,0]
complist[4] = [-1,-1,1]




for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

    criteria = np.where((Stars['R']>(binprefacs[irad]*minrad)) & (Stars['R']<(binprefacs[irad]*maxrad)))[0]

    for cnum in [0,1,2]:


        for ikey,key in enumerate(pltkeys):
            if 'inv' in key:
                median = np.sqrt(1./np.nanmedian(CStats[cnum][key]))
                lo = np.sqrt(1./np.nanpercentile(CStats[cnum][key],14.)) - median
                hi = np.sqrt(1./np.nanpercentile(CStats[cnum][key],86.)) - median

            elif key=='alpha':
                median = (180./np.pi)*np.nanmedian(CStats[cnum][key])
                lo = (180./np.pi)*np.nanpercentile(CStats[cnum][key],14.) - median
                hi = (180./np.pi)*np.nanpercentile(CStats[cnum][key],86.) - median
            elif key=='f':
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median

            else: # Lx, Ly, Lz
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median

            radval = np.nanmedian(Stars['R'][criteria])
            if comptag[minrad][cnum]=='bar':
                axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='black')
                axlist[ikey].scatter(radval,median,marker='o',edgecolor='none',facecolor='black')
            elif comptag[minrad][cnum]=='disc':
                axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='blue')
                axlist[ikey].scatter(radval,median,marker='o',edgecolor='none',facecolor='blue')
            elif comptag[minrad][cnum]=='knot':
                axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='red')
                axlist[ikey].scatter(radval,median,marker='o',edgecolor='none',facecolor='red')






for ax in axlist:
    ax.tick_params(axis="both",direction="in",which="both")

ax1.axis([0.,5.,0.,1.1])
ax2.axis([0.,5.,-50,50.])
ax3.axis([0.,5.,-50,50.])

# add the Lz curves
ax4.plot(modelradii,-np.array(barlz),color='black')
ax4.plot(modelradii,-np.array(potentialmeasuredlz),color='blue')
ax4.plot(modelradii,-np.array(massenclosedlz),color='blue',linestyle='dashed')

ax4.axis([0.,5.,-1000,100.])


ax5.axis([0.,5.,0.,90.])
ax6.axis([0.,5.,0.,300.])
ax7.axis([0.,5.,0.,300.])
ax8.axis([0.,5.,0.,300.])

for ax in [ax1,ax2,ax3,ax4]:
    ax.set_xticklabels(())

ax6.set_xlabel('radius (kpc)',x=1.3)

for ikey,key in enumerate(pltkeys):
    axlist[ikey].set_ylabel(key)

plt.savefig('figures/fitvalues_{}.png'.format(modeltag),dpi=300)
