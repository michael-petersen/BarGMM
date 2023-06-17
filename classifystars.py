
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
#from models.bulgemock import *
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


print('{0:4s}{1:3s}{2:9s}{3:9s} {4:9s} {5:9s} {6:9s} {7:9s} {8:9s} {9:9s}'.format(' ',' ','f','Lx','Ly','Lz','alpha','sx','sy','sz'))

# markdown version for GitHub viewing

f = open(inputdir+'fits/README.md','w')
print('|comp|radii| f | L<sub>x</sub> | L<sub>y</sub> | L<sub>z</sub> | angle | w<sub>x</sub> | w<sub>y</sub> | w<sub>z</sub> |',file=f)
print('|---|---|---| ---| --- | ---| --- | --- | --- | --- |',file=f)


g = open(inputdir+'fits/table.tex','w')
print('comp &radii & f & $L_x$& $L_y$ & $L_z$ & $\alpha$ & $\sigma_x$ & $\sigma_y$ & $\sigma_z$ \\\\',file=g)


h = open(inputdir+'fits/table.csv','w')
print('comp,binmin,binmax,starmed,starmed-,starmed+,f,f-,f+,Lx,Lx-,Lx+,Ly,Ly-,Ly+,Lz,Lz-,Lz+,alpha,alpha-,alpha+,sigmax,sigmax-,sigmax+,sigmay,sigmay-,sigmay+,sigmaz,sigmaz-,sigmaz+,',file=h)

#for irad,rads in enumerate():
for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}".format(modeltag,minrad,maxrad)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)
    print('{0:2.1f}-{1:2.1f}'.format(np.round(binprefacs[irad]*minrad,1),np.round(binprefacs[irad]*maxrad,1)))

    criteria = np.where((Stars['R']>(binprefacs[irad]*minrad)) & (Stars['R']<(binprefacs[irad]*maxrad)))[0]

    for cnum in [0,1,2]:

        # print to markdown
        print('|'+comptag[minrad][cnum],end='|',file=f)
        print('{0:2.1f}-{1:2.1f}'.format(np.round(binprefacs[irad]*minrad,1),np.round(binprefacs[irad]*maxrad,1)),end='|',file=f)

        # print to latex
        print(comptag[minrad][cnum],end='&',file=g)
        print('{0:2.1f}-{1:2.1f}'.format(np.round(binprefacs[irad]*minrad,1),np.round(binprefacs[irad]*maxrad,1)),end='&',file=g)

        # print to csv
        print(comptag[minrad][cnum],end=',',file=h)
        print('{0:2.1f},{1:2.1f}'.format(np.round(binprefacs[irad]*minrad,1),np.round(binprefacs[irad]*maxrad,1)),end=',',file=h)
        print('{0:3.2f},{1:3.2f},{2:3.2f}'.format(np.nanmedian(Stars['R'][criteria]),np.nanpercentile(Stars['R'][criteria],14.),np.nanpercentile(Stars['R'][criteria],86.)),end=',',file=h)

        for ikey,key in enumerate(pltkeys):
            if 'inv' in key:
                median = np.sqrt(1./np.nanmedian(CStats[cnum][key]))
                lo = np.sqrt(1./np.nanpercentile(CStats[cnum][key],14.)) - median
                hi = np.sqrt(1./np.nanpercentile(CStats[cnum][key],86.)) - median
                print(' {0:9.4f}'.format(median),end='')
                print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='|',file=f)
                print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='&',file=g)
                print('{0},-{2},+{1}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end=',',file=h)
            elif key=='alpha':
                median = (180./np.pi)*np.nanmedian(CStats[cnum][key])
                lo = (180./np.pi)*np.nanpercentile(CStats[cnum][key],14.) - median
                hi = (180./np.pi)*np.nanpercentile(CStats[cnum][key],86.) - median
                print(' {0:9.4f}'.format(median),end='')
                print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='|',file=f)
                print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='&',file=g)
                print('{0},-{2},+{1}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end=',',file=h)
            elif key=='f':
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median
                print(' {0:9.4f}'.format(median),end='')
                print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,4),np.round(np.abs(lo),4),np.round(np.abs(hi),4)),end='|',file=f)
                print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,4),np.round(np.abs(lo),4),np.round(np.abs(hi),4)),end='&',file=g)
                print('{0},-{2},+{1}'.format(np.round(median,4),np.round(np.abs(lo),4),np.round(np.abs(hi),4)),end=',',file=h)
            else: # Lx, Ly, Lz
                median = np.nanmedian(CStats[cnum][key])
                lo = np.nanpercentile(CStats[cnum][key],14.) - median
                hi = np.nanpercentile(CStats[cnum][key],86.) - median
                print(' {0:9.4f}'.format(median),end='')
                print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='|',file=f)
                print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end='&',file=g)
                print('{0},-{2},+{1}'.format(np.round(median,1),np.round(np.abs(lo),1),np.round(np.abs(hi),1)),end=',',file=h)

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
        print(' '+comptag[minrad][cnum])
        print('',file=f)
        print('\\\\',file=g)
        print('',file=h)

f.close()

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

plt.savefig('figures/fitvalues_{}.png'.format(modeltag),dpi=300)

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
            f = h5py.File(inputdir+"classifications/AllClassifications_{}.h5".format(modeltag),"w")

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
