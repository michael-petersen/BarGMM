
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
from src.tableprint import *

#from models.apogee import *
#from models.bulgemock import *
from models.barmock import *

Stars = read_mock_file(datafile)


classify = False


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

A = table_print(format='terminal')
B = table_print(format='markdown',outputfile=open(inputdir+'fits/README{0}.md'.format(appendix),'w'))
C = table_print(format='csv',outputfile=open(inputdir+'fits/table{0}.csv'.format(appendix),'w'))
D = table_print(format='tex',outputfile=open(inputdir+'fits/table{0}.tex'.format(appendix),'w'))

#for irad,rads in enumerate():
for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

    criteria = np.where((Stars['R']>(binprefacs[irad]*minrad)) & (Stars['R']<(binprefacs[irad]*maxrad)))[0]

    A.print_column1(comptag,binprefacs,minrad,maxrad,irad,0,Stars,criteria)

    for cnum in [0,1,2]:

        for p in [B,C,D]:
            p.print_column1(comptag,binprefacs,minrad,maxrad,irad,cnum,Stars,criteria)

        for ikey,key in enumerate(pltkeys):
            if 'inv' in key:
                median = np.sqrt(1./np.nanmedian(CStats[cnum][key]))
                lo = np.sqrt(1./np.nanpercentile(CStats[cnum][key],14.)) - median
                hi = np.sqrt(1./np.nanpercentile(CStats[cnum][key],86.)) - median
                for p in [A,B,C,D]:
                    p.print_key_column(median,lo,hi,rounding=1)

            elif key=='alpha':
                median = (180./np.pi)*np.nanmedian(CStats[cnum][key])
                lo = (180./np.pi)*np.nanpercentile(CStats[cnum][key],14.) - median
                hi = (180./np.pi)*np.nanpercentile(CStats[cnum][key],86.) - median
                for p in [A,B,C,D]:
                    p.print_key_column(median,lo,hi,rounding=1)
            elif key=='f':
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median
                for p in [A,B,C,D]:
                    p.print_key_column(median,lo,hi,rounding=4)
            else: # Lx, Ly, Lz
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median
                for p in [A,B,C,D]:
                    p.print_key_column(median,lo,hi,rounding=1)

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

        for p in [A,B,C,D]:
            p.print_columnX(comptag,minrad,cnum)

for p in [B,C,D]:
    p.f.close()

for ax in axlist:
    ax.tick_params(axis="both",direction="in",which="both")

ax1.axis([0.,5.,0.,1.1])
ax2.axis([0.,5.,-50,50.])
ax3.axis([0.,5.,-50,50.])
ax4.axis([0.,5.,-800,100.])
ax5.axis([0.,5.,0.,90.])
ax6.axis([0.,5.,0.,220.])
ax7.axis([0.,5.,0.,220.])
ax8.axis([0.,5.,0.,220.])

for ax in [ax1,ax2,ax3,ax4]:
    ax.set_xticklabels(())

ax6.set_xlabel('radius (kpc)',x=1.3)

for ikey,key in enumerate(pltkeys):
    axlist[ikey].set_ylabel(key)

plt.savefig('figures/fitvalues_{0}{1}.png'.format(modeltag,appendix),dpi=300)

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
        #criteria = np.where((Stars['R']>(minrad)) & (Stars['R']<(maxrad)))[0]
        criteria = np.where((Stars['R']>(binprefacs[irad]*minrad)) & (Stars['R']<(binprefacs[irad]*maxrad)))[0]

        # do the probabilistic classification
        # allprobs is the full set of classifications
        # percentileprob is the 50h percentile
        # errorprob is the 1-sigma error on the classification
        allprobs,percentileprob,errorprob = make_all_probabilities(Stars,criteria,CStats,nchains=20)


        # which component is which?
        disccomp,barcomp,knotcomp = compnum[minrad]

        if minrad==radii[0][0]:
            f = h5py.File(inputdir+"classifications/AllClassifications_{0}{1}.h5".format(modeltag,appendix),"w")

        nsamples = 20
        for indx,starnum in enumerate(criteria):
            probabilities = np.zeros([nsamples,4]) # always disc, bar, knot, [x,y,z,Lx,Ly,Lz]
            if disccomp>=0: probabilities[:,0] = allprobs[:,indx,disccomp]
            if barcomp>=0:  probabilities[:,1] = allprobs[:,indx,barcomp]
            if knotcomp>=0: probabilities[:,2] = allprobs[:,indx,knotcomp]
            probabilities[0:7,3] = [Stars['R'][starnum],Stars['x'][starnum],Stars['y'][starnum],Stars['z'][starnum],Stars['Lx'][starnum],Stars['Ly'][starnum],Stars['Lz'][starnum]]
            if binprefacs[irad] > 0.5:
                try:
                    dset = f.create_dataset(Stars['apogee_id'][starnum], data=probabilities)
                except:
                    print(Stars['apogee_id'][starnum])
            else:
                dset = f.create_dataset(Stars['apogee_id'][starnum]+b'*', data=probabilities)

        if minrad==radii[-1][0]: f.close()



        # print the classification
        #radii = [minrad,maxrad]
        #printdir = inputdir+'classifications/'
        #print_classification(Stars,criteria,radii,percentileprob,errorprob,disccomp,barcomp,knotcomp,printdir,mockanalysis)
