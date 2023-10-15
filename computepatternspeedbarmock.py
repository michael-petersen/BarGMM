import matplotlib.pyplot as plt
# %matplotlib inline
import numpy as np
from astropy.io import fits
import tqdm
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
rc('text')#, usetex=True)

import nmmn.plots
import vorbin

wolfram=nmmn.plots.wolframcmap() # for Mathematica's cmap
parula=nmmn.plots.parulacmap() # for MATLAB's cmap
turbo=nmmn.plots.turbocmap() # Turbo
import warnings
warnings.filterwarnings('ignore')

import pandas as pd

import h5py

from scipy.optimize import curve_fit

from scipy.optimize import fmin, minimize


#appendix="_apog0"
appendix="_apog1"

# load the classification data
f = h5py.File("data/barmock/classifications/AllClassifications_ModelBarYoungMock10000b{}.h5".format(appendix),"r")
median_knot_probability = []
median_bar_probability = []
median_disc_probability = []
std_knot_probability = []
std_bar_probability = []
std_disc_probability = []
all_knot_probability = []
all_bar_probability = []
all_disc_probability = []
RxyzL = []


for key in tqdm.tqdm_notebook(f.keys()):
    median_knot_probability.append(np.nanmedian(f[key][:,2]))
    median_bar_probability.append(np.nanmedian(f[key][:,1]))
    median_disc_probability.append(np.nanmedian(f[key][:,0]))
    std_knot_probability.append(np.nanstd(f[key][:,2]))
    std_bar_probability.append(np.nanstd(f[key][:,1]))
    std_disc_probability.append(np.nanstd(f[key][:,0]))
    all_knot_probability.append((f[key][:,2]))
    all_bar_probability.append((f[key][:,1]))
    all_disc_probability.append((f[key][:,0]))
    RxyzL.append(f[key][:,3][:12])
f.close() # don't forget to close the file!

mask_knot2 = (np.array(median_knot_probability)>0.8)&(np.array(std_knot_probability)<0.15)
mask_bar2 = (np.array(median_bar_probability)>0.8)&(np.array(std_bar_probability)<0.15)
mask_disc2 = (np.array(median_disc_probability)>0.8)&(np.array(std_disc_probability)<0.15)
RxyzL_ = np.array(RxyzL)
R = RxyzL_[:,0]
x = RxyzL_[:,1]
y = RxyzL_[:,2]
z = RxyzL_[:,3]
Lx = RxyzL_[:,7]
Ly = RxyzL_[:,8]
Lz = RxyzL_[:,9]*-1
eR = RxyzL_[:,10]
eLz = RxyzL_[:,11]

# find the median bar radius per bin per classification
radvals = np.zeros(4) # we know there are four independent bins with bar classifications
criteria = (R>0) & (R<1)
radvals[0] = np.nansum(R[criteria]*np.array(median_bar_probability)[criteria])/np.nansum(np.array(median_bar_probability)[criteria])

criteria = (R>1) & (R<2)
radvals[1] = np.nansum(R[criteria]*np.array(median_bar_probability)[criteria])/np.nansum(np.array(median_bar_probability)[criteria])

criteria = (R>2) & (R<3)
radvals[2] = np.nansum(R[criteria]*np.array(median_bar_probability)[criteria])/np.nansum(np.array(median_bar_probability)[criteria])

criteria = (R>3) & (R<4)
radvals[3] = np.nansum(R[criteria]*np.array(median_bar_probability)[criteria])/np.nansum(np.array(median_bar_probability)[criteria])


# load the csv data, but only the independent bins!
data = pd.read_csv('data/barmock/fits/table{}.csv'.format(appendix))

mask_disc = np.where(data['comp']=='disc')
mask_knot = np.where(data['comp']=='knot')
mask_bar = np.where(data['comp']=='bar')

mask_bar = [np.array([1,4,7,11])]


# we're just fitting a linear model, no intercept.
def func(x,a):
    return a*x

popt, pcov = curve_fit(func, (radvals)**2, data['Lz'].values[mask_bar],sigma=data['sigmaz'].values[mask_bar])

print('The pattern speed is {0}pm{1}'.format(popt,np.sqrt(pcov[0])))



# this is the all-probability-fit routine:
# this method is subject to a lot of measurement uncertainty at 20 samples
# if we wanted to use this method, we would need to improve the sampling
weights = np.array(all_bar_probability)

def linear(xs,m):
    return m*xs

res = []
for indx, i in tqdm.tqdm_notebook(enumerate(weights.T)):
    def loss2(w):
        XS,YS = (R+np.random.normal(0.,0.5*eR))**2,Lz+np.random.normal(0.,0.5*eLz) # Lz uncertainties too correlated in the mock, reduce overall uncertainties as a test to compensate
        popt, pcov = curve_fit(linear,XS,YS)
        ypred = linear(XS,w)
        # here we will use a weighted least squares
        return np.sum(i*(YS-ypred)**2)

    res.append(minimize(loss2, popt, method='nelder-mead'))


slopes = []
for indx, i in enumerate(res):
    slopes.append(i['x'])
slopes = np.array(slopes).T[0]

print(np.mean(slopes))
print(np.std(slopes))
