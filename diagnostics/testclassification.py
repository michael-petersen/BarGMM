
# basic imports
import numpy as np
from numpy.linalg import eig, inv

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


from src.localtools import *
from src.fitmanagement import *

from models.apogee import *

Stars = read_mock_file(datafile)



minrad,maxrad = 15,25
binprefac = 0.1

minrad,maxrad = 1,2
binprefac = 1.

# open the correct chain
directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}".format(modeltag,minrad,maxrad)
inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

# define the stars we want to classify
criteria = np.where((Stars['R']>(binprefac*minrad)) & (Stars['R']<(binprefac*maxrad)))[0]





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


nsamples = 20
for indx,starnum in enumerate(criteria):
    probabilities = np.zeros([nsamples,4]) # always disc, bar, knot, [x,y,z,Lx,Ly,Lz]
    if disccomp>=0: probabilities[:,0] = allprobs[:,indx,disccomp]
    if barcomp>=0:  probabilities[:,1] = allprobs[:,indx,barcomp]
    if knotcomp>=0: probabilities[:,2] = allprobs[:,indx,knotcomp]
    probabilities[0:7,3] = [Stars['R'][starnum],Stars['x'][starnum],Stars['y'][starnum],Stars['z'][starnum],Stars['Lx'][starnum],Stars['Ly'][starnum],Stars['Lz'][starnum]]
    dset = f.create_dataset(Stars['apogee_id'][starnum], data=probabilities)


f.close()


f = h5py.File("tst.h5","r")
#f.keys() # these are the APOGEE IDs

for key in f.keys():
    if f[key][0,1] > 0.9:
        #print(f[key][6,3])
        _ = plt.scatter(f[key][1,3],f[key][2,3],color='black',s=1.)



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
