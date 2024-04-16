
# basic imports
import numpy as np
from numpy.linalg import eig, inv
import h5py

# plotting elements
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
mpl.rcParams['font.weight'] = 'medium';mpl.rcParams['xtick.labelsize'] = 10;mpl.rcParams['ytick.labelsize'] = 10
plt.ion()

import scipy.interpolate as interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import scipy

#import sys; sys.path.append("../")

from src.localtools import *
from src.fitmanagement import *

from models.apogee import *
from models.barmock import *

Stars = read_mock_file(datafile)



minrad,maxrad = 15,25
binprefac = 0.1

minrad,maxrad = 1,2
binprefac = 1.

disccomp,barcomp,knotcomp = compnum[minrad]

# open the correct chain
directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

# define the stars we want to classify
criteria = np.where((Stars['R']>(binprefac*minrad)) & (Stars['R']<(binprefac*maxrad)) & (Stars['apogee']==1))[0]

allprobs,percentileprob,errorprob = make_all_probabilities(Stars,criteria,CStats,nchains=10)

# let's check the bar membership against the classifications
# allprobs.shape
# (10, 2643, 3)

# this size works well for single-column journal figures
fig = plt.figure(figsize=(3.87,3.0),facecolor='white')

xmin = 0.17
ymin = 0.13
dx = 0.65
dy = 0.83

classedsample = Stars['x'][criteria][(Stars['bulge'][criteria]>=0)].size # the size of the total sample

# indices of the classed sample, could look at their actual membership logs
Stars['apogee_id'][criteria][(Stars['bulge'][criteria]>=0)]

barmems = ((np.nanmedian(allprobs[:,:,barcomp],axis=0)<0.5) & (Stars['bulge'][criteria]>0.)& (Stars['bulge'][criteria]>=0))
print('particles that are a part of the bar but not classed as such (WRONG)')
print(Stars['x'][criteria][barmems].size) # particles that are a part of the bar but not classed as such (WRONG)

barmems = ((np.nanmedian(allprobs[:,:,barcomp],axis=0)<0.5) & (Stars['bulge'][criteria]<=0.)& (Stars['bulge'][criteria]>=0))
print(' particles that are not a part of the bar and are not classed as such (CORRECT)')
print(Stars['x'][criteria][barmems].size) # particles that are not a part of the bar and are not classed as such (CORRECT)

barmems = ((np.nanmedian(allprobs[:,:,barcomp],axis=0)>0.5) & (Stars['bulge'][criteria]<=0.)& (Stars['bulge'][criteria]>=0))
print('particles that are not a part of the bar and are classed as such (WRONG)')
print(Stars['x'][criteria][barmems].size) # particles that are not a part of the bar and are classed as such (WRONG)

barmems = ((np.nanmedian(allprobs[:,:,barcomp],axis=0)>0.5) & (Stars['bulge'][criteria]>0.)& (Stars['bulge'][criteria]>=0))
print('particles that are a part of the bar and are classed as such (CORRECT)')
print(Stars['x'][criteria][barmems].size) # particles that are a part of the bar and are classed as such (CORRECT)
# disc stars are rarely classed as bar stars (3% contamination)

barmems = ((np.nanmedian(allprobs[:,:,disccomp],axis=0)<0.5) & (Stars['bulge'][criteria]>0.)& (Stars['bulge'][criteria]>=0))
print(Stars['x'][criteria][barmems].size) # particles that are a part of the bar but not classed as such (CORRECT)
barmems = ((np.nanmedian(allprobs[:,:,disccomp],axis=0)<0.5) & (Stars['bulge'][criteria]<=0.)& (Stars['bulge'][criteria]>=0))
print(Stars['x'][criteria][barmems].size) # particles that are not a part of the disc and are not classed as such (WRONG)
barmems = ((np.nanmedian(allprobs[:,:,disccomp],axis=0)>0.5) & (Stars['bulge'][criteria]<=0.)& (Stars['bulge'][criteria]>=0))
print(Stars['x'][criteria][barmems].size) # particles that are not a part of the bar and are classed as such (CORRECT)
barmems = ((np.nanmedian(allprobs[:,:,disccomp],axis=0)>0.5) & (Stars['bulge'][criteria]>19.)& (Stars['bulge'][criteria]>=0))
print(Stars['x'][criteria][barmems].size) # particles that are a part of the bar and are classed as disc (WRONG)
# the most common failure mode is bar-classified stars leaking into the disc component. these tend to be the larger Lz stars, where the classification is more ambigious. However, these stars are not enough to bias the results in the mock dataset, which recovers the correct pattern speed (derived from the mean Lz of the bar stars).

barmems = ((np.nanmedian(allprobs[:,:,knotcomp],axis=0)>0.5) & (Stars['bulge'][criteria]<10.)& (Stars['bulge'][criteria]>=0))
print(Stars['x'][criteria][barmems].size) # particles that are not a part of the bar and are classed as knot (WRONG)
barmems = ((np.nanmedian(allprobs[:,:,knotcomp],axis=0)>0.5) & (Stars['bulge'][criteria]>10.)& (Stars['bulge'][criteria]>=0))
print(Stars['x'][criteria][barmems].size) # particles that are a part of the bar and are classed as knot (CORRECT)
# knot stars are classed as bar in PWK21.



ax1 = fig.add_axes([xmin+0*dx     ,ymin+0*dy,dx,dy])   # main figure
ax2 = fig.add_axes([xmin+1*dx+0.02,ymin+0*dy,0.02,dy]) # colourbar


barmems = ((Stars['bulge'][criteria]>10.)& (Stars['bulge'][criteria]>=0))
ax1.scatter(np.nanmedian(allprobs[:,barmems,barcomp],axis=0),Stars['Lz'][criteria][barmems],facecolor='firebrick',edgecolor='None',s=3.)

ax1.scatter(np.nanmedian(allprobs[:,barmems,disccomp],axis=0),Stars['Lz'][criteria][barmems],facecolor='limegreen',edgecolor='None',s=3.)



ax1.scatter(Stars['x'][criteria][barmems],Stars['y'][criteria][barmems],facecolor='firebrick',edgecolor='None',s=3.)

discmems = ((np.nanmedian(allprobs[:,:,disccomp],axis=0)>0.7) & (Stars['bulge'][criteria]>19.))
#ax1.scatter(Stars['x'][criteria][discmems],Stars['y'][criteria][discmems],facecolor='limegreen',edgecolor='None',s=3.)

ax1.axis([-10.,10.,-10.,10.])
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

ax1.tick_params(axis="both",direction="in",which="both")

cmin,cmax = 0.,1.
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
cmap = cm.magma
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,norm=norm)
cb1.set_label('colourbar label')
cb1.set_ticks([0.,0.5,1.0])
cb1.ax.minorticks_off()

plt.savefig('/Users/mpetersen/Downloads/comp_bar.png',dpi=300)




ax1.scatter(np.nanmedian(allprobs[:,:,barcomp],axis=0),Stars['bulge'][criteria],facecolor='firebrick',edgecolor='None',s=3.)
ax1.scatter(np.nanmedian(allprobs[:,:,disccomp],axis=0),Stars['bulge'][criteria],facecolor='limegreen',edgecolor='None',s=3.)
ax1.scatter(np.nanmedian(allprobs[:,:,knotcomp],axis=0),Stars['bulge'][criteria],facecolor='indigo',edgecolor='None',s=3.)
ax1.axis([0.,1.,-1.,21.])
ax1.set_xlabel('median bar probability')
ax1.set_ylabel('bar membership')




for ii in range(9,10):
    indx = criteria[ii]
    L = np.cov(np.array([Stars['eLx'][indx],Stars['eLy'][indx],Stars['eLz'][indx]]))
    print(L)
    for cnum in [0,1,2]:
        cosa,sina = np.cos(COMPS[cnum][7]),np.sin(COMPS[cnum][7])
        S = np.array([[COMPS[cnum][4]*cosa*cosa + COMPS[cnum][5]*sina*sina,(COMPS[cnum][4]-COMPS[cnum][5])*cosa*sina,0],\
                         [(COMPS[cnum][4]-COMPS[cnum][5])*cosa*sina,COMPS[cnum][4]*sina*sina + COMPS[cnum][5]*cosa*cosa,0],\
                         [0,0,COMPS[cnum][6]]])
        C = L + S
        print(indx,cnum,np.linalg.det(C),det_by_hand(C[0][0],C[0][1],C[0][2],C[1][1],C[1][2],C[2][2]))


def det_by_hand(a,b,c,e,f,i):
    return a*e*i - a*f*f - e*c*c - i*b*b + 2*b*c*f


indx = 1756
for key in Stars.keys():
    print(key,Stars[key][indx])



# do the probabilistic classification
# allprobs is the full set of classifications
# percentileprob is the 50h percentile
# errorprob is the 1-sigma error on the classification
allprobs,percentileprob,errorprob = make_all_probabilities(Stars,criteria,CStats,nchains=10,nonorm=False)
allprobs2,percentileprob2,errorprob2 = make_all_probabilities(Stars,criteria,CStats,nchains=10,nonorm=True)

# allprobs.shape: [nchains,criteria.size,3] where the 3 is

minrad,maxrad = radii[0],radii[1]

Stars['discprob'] = np.zeros(Stars['x'].size)
Stars['ediscprob'] = np.zeros(Stars['x'].size)

if disccomp>=0:
    Stars['discprob'][criteria] = percentileprob[:,disccomp]
    Stars['ediscprob'][criteria] = errorprob[:,disccomp]


Stars['barprob'] = np.zeros(Stars['x'].size)
Stars['ebarprob'] = np.zeros(Stars['x'].size)

if barcomp>=0:
    Stars['barprob'][criteria] = percentileprob[:,barcomp]
    Stars['ebarprob'][criteria] = errorprob[:,barcomp]


Stars['knotprob'] = np.zeros(Stars['x'].size)
Stars['eknotprob'] = np.zeros(Stars['x'].size)

if knotcomp>=0:
    Stars['knotprob'][criteria] = percentileprob[:,knotcomp]
    Stars['eknotprob'][criteria] = errorprob[:,knotcomp]



import h5py

disccomp,barcomp,knotcomp = compnum[minrad]


f = h5py.File("data/apogee/classifications/AllClassifications_APOGEE_all_feh8_SNRc.h5","r")

stride = 10
stridecounter = 0
for key in f.keys():
    if (stridecounter % stride ) == 0:
        #if f[key][0,1] > 0.8:
        #    _ = plt.scatter(f[key][1,3],f[key][2,3],color='black',s=1.)
        if f[key][0,0] > 0.8:
            _ = plt.scatter(f[key][1,3],f[key][2,3],color='blue',s=1.)
    stridecounter += 1

barr = np.linspace(-3.5,3.5,100)
plt.plot(barr*np.cos(-10*np.pi/180.),barr*np.sin(-10*np.pi/180.),color='black',lw=2.)
plt.plot(barr*np.cos(-20*np.pi/180.),barr*np.sin(-20*np.pi/180.),color='black',lw=2.)
plt.plot(barr*np.cos(-30*np.pi/180.),barr*np.sin(-30*np.pi/180.),color='black',lw=2.)


# need a way to match these against the apogee ids


f = h5py.File("data/bulgemock/classifications/AllClassifications_ModelBarYoungBulgeMock10000.h5","r")
g = h5py.File("data/bulgemock/classifications/AllClassifications_ModelBarYoungBulgeMock10000_apog0.h5","r")

# match up by index
keyarray = np.array([int(key.split('.')[0]) for key in f.keys()])
xy,xind,yind = np.intersect1d(keyarray,Stars['index'].astype('int'),return_indices=True)

keymax = 1000
for ikey,key in enumerate(f.keys()):
    if ikey<keymax:
        #if f[key][0,3]<1.0:
        if Stars['apogee'][yind][ikey]==1:
            #print(f[key][6,3])
            _ = plt.scatter(f[key][:,1],g[key][:,1],color='black',s=1.)


plt.xlabel('apogee==1')
plt.ylabel('apogee==0')





f = h5py.File("/Users/mpetersen/Notebooks/Data/BarGMM/data/apogee/classifications/AllClassifications.h5","r")

median_bar_probability = []
for key in f.keys():
    median_bar_probability.append(np.nanmedian(f[key][:,1]))



xbar,ybar = [],[]
for key in f.keys():
    if f[key][0,0] > 0.9:
        #print(f[key][6,3])
        xbar.append(f[key][1,3])
        ybar.append(f[key][2,3])


# Stars.keys()
# dict_keys(['l', 'b', 'dist', 'edist', 'pml', 'epml', 'pmb', 'epmb', 'pmlpmbcorr', 'vlos', 'evlos', 'apogee', 'bulge', 'feh', 'alpha', 'age', 'index', 'apogee_id', 'u', 'v', 'w', 'x', 'y', 'z', 'r', 'R', 'Lx', 'Ly', 'Lz', 'L', 'eLx', 'eLy', 'eLz'])

# find bulge stars with very large Ly
#

# Ly vs probability of bulge
plt.scatter(percentileprob[:,2],Stars['Ly'][criteria],color='black',s=4.)



membership = (percentileprob[:,2]>0.8)
plt.scatter(Stars['Lz'][criteria],Stars['Ly'][criteria],color='grey',s=2.)
plt.scatter(Stars['Lz'][criteria][membership],Stars['Ly'][criteria][membership],color='black',s=4.)



for m in range(0,len(Stars['Lz'][criteria][membership])):
    x,y = Stars['Lz'][criteria][membership][m],Stars['Ly'][criteria][membership][m]
    dx,dy = np.nanstd(Stars['eLz'][criteria][membership][m]),np.nanstd(Stars['eLy'][criteria][membership][m])
    _ = plt.plot([x,x],[y-dy,y+dy],color='black',lw=0.5)
    _ = plt.plot([x-dx,x+dx],[y,y],color='black',lw=0.5)


membership = (percentileprob[:,2]>0.9)

weirdones = np.where(np.abs(Stars['Ly'][criteria][membership])>100)[0]
plt.scatter(Stars['Lz'][criteria][membership][weirdones],Stars['Ly'][criteria][membership][weirdones],color='red',s=4.)

print(Stars['Ly'][criteria][membership][weirdones])
print(percentileprob[:,2][membership][weirdones])


tstmatrix = np.array([[1,2,3],[4,5,6],[7,8,9]])
tstmatrix = np.random.rand(3,3,dtype=np.float64)
C = np.cov(tstmatrix)
Cn = np.corrcoef(tstmatrix)

ex,ey,ez = np.diag(C) # get the individual data values
exd = ex
eyd = ey
ezd = ez
ex,ey,ez = np.sqrt(ex),np.sqrt(ey),np.sqrt(ez)

Cf = np.zeros_like(Cn)
Cf[0][0],Cf[0][1],Cf[0][2] = exd   ,ex*ey*Cn[0][1],ex*ez*Cn[0][2]
Cf[1][0],Cf[1][1],Cf[1][2] = ey*ex*Cn[1][0],eyd   ,ey*ez*Cn[1][2]
Cf[2][0],Cf[2][1],Cf[2][2] = ez*ex*Cn[2][0],ez*ey*Cn[2][1],ezd

print(Cf)
print(C)

# now try determinant of C by hand
print(np.linalg.det(C))

def det_by_hand(a,b,c,e,f,i):
    return a*e*i - a*f*f - e*c*c - i*b*b + 2*b*c*f

print(det_by_hand(C[0][0],C[0][1],C[0][2],C[1][1],C[1][2],C[2][2]))
