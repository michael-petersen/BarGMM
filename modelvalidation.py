

# basic imports
import numpy as np

# data storage imports
import h5py

# local tools
from src.localtools import *
from src.fitmanagement import *
from src.tableprint import *
from src.parameterfigure import *

# select the model of interest
from models.barmock import *

import pandas as pd



plt.figure(facecolor='white')

# if we were to assume this was all circular velocity...
plt.plot(modelradii**2,potentialmeasuredlz,color='blue')
plt.plot(modelradii**2,massenclosedlz,color='blue',linestyle='dashed')
plt.plot(modelradii**2,220*modelradii,color='blue',linestyle='dotted')

width1 = 'Lz-'
width2 = 'Lz+'

width1 = 'sigmaz'
width2 = 'sigmaz'

# add the fits
F = pd.read_csv('data/barmock/fits/table_apog1.csv')
#F = pd.read_csv('data/bulgemock/fits/table_apog1.csv')
discselect = F.loc[F['comp']=='disc']
plt.scatter(discselect['starmed']**2,np.abs(discselect['Lz']),edgecolor='blue',facecolor='blue',s=30.)
for i in range(0,discselect['starmed'].values.size):
    plt.plot([discselect['starmed'].values[i]**2,discselect['starmed'].values[i]**2],\
             [np.abs(discselect['Lz']-np.abs(discselect[width1])).values[i],np.abs(discselect['Lz']+np.abs(discselect[width2])).values[i]],color='blue')


G = pd.read_csv('data/barmock/fits/table_apog0.csv')
#G = pd.read_csv('data/bulgemock/fits/table_apog1.csv')
discselect = G.loc[G['comp']=='disc']
plt.scatter(discselect['starmed']**2,np.abs(discselect['Lz']),edgecolor='blue',facecolor='white',s=30.)
for i in range(0,discselect['starmed'].values.size):
    plt.plot([discselect['starmed'].values[i]**2,discselect['starmed'].values[i]**2],\
             [np.abs(discselect['Lz']-np.abs(discselect[width1])).values[i],np.abs(discselect['Lz']+np.abs(discselect[width2])).values[i]],color='blue',linestyle='dotted',lw=3.)



# the mass enclosed is an overestimate of the angular momentum at all radii
# there can be scatter around the relation from stars oscillating around their guiding centre
# there is always a bias below circular orbits, since orbits have radial energy

plt.plot(modelradii**2,barlz,color='black')

barselect = F.loc[F['comp']=='bar']
plt.scatter(barselect['starmed']**2,np.abs(barselect['Lz']),edgecolor='black',facecolor='black',s=30.)
for i in range(0,barselect['starmed'].values.size):
    plt.plot([barselect['starmed'].values[i]**2,barselect['starmed'].values[i]**2],\
             [np.abs(barselect['Lz']-np.abs(barselect[width1])).values[i],np.abs(barselect['Lz']+np.abs(barselect[width2])).values[i]],color='black')


barselect = G.loc[G['comp']=='bar']
plt.scatter(barselect['starmed']**2,np.abs(barselect['Lz']),edgecolor='black',facecolor='white',s=30.)
for i in range(0,barselect['starmed'].values.size):
    plt.plot([barselect['starmed'].values[i]**2,barselect['starmed'].values[i]**2],\
             [np.abs(barselect['Lz']-np.abs(barselect[width1])).values[i],np.abs(barselect['Lz']+np.abs(barselect[width2])).values[i]],color='black',linestyle='dotted',lw=3.)




plt.axis([0.,85.,0.,2500.])
#plt.axis([0.,18.,0.,1200.])

plt.xlabel('square radius (kpc$^2$)')
plt.ylabel('angular momentum (kpc km/s)')

plt.tight_layout()
plt.savefig('figures/angularmomentum_diagnostic.png',dpi=300)
