

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

#plt.plot(np.linspace(0.,5.,100),24*(np.linspace(0.,5.,100)**2),color='black')
#plt.plot(np.linspace(0.,5.,100),33*(np.linspace(0.,5.,100)**2),color='black',linestyle='dotted')

medval = 'Lz'
width1 = 'Lz-'
width2 = 'Lz+'

width1 = 'sigmaz'
width2 = 'sigmaz'


medval = 'alpha'
width1 = 'alpha-'
width2 = 'alpha+'

"""

medval = 'f'
width1 = 'f-'
width2 = 'f+'

medval = 'sigmaz'
width1 = 'sigmaz-'
width2 = 'sigmaz+'
"""

# add the fits
F = pd.read_csv('data/apogee/fits/table.csv')
G = pd.read_csv('data/apogee/fits/table_astroNN.csv')
H = pd.read_csv('data/apogee/fits/table_StarHorse.csv')

modsymbols = ['*','o','^']
# Bailer-Jones are stars
# AstroNN are circles
# StarHorse are triangles
plotdisc = True
plotbar = True
plotknot = False

for imod,mod in enumerate([F,G,H]):

    if plotdisc:
        discselect = mod.loc[mod['comp']=='disc']
        plt.scatter(discselect['starmed'],np.abs(discselect[medval]),marker=modsymbols[imod],edgecolor='blue',facecolor='blue',s=30.)
        for i in range(0,discselect['starmed'].values.size):
            plt.plot([discselect['starmed'].values[i],discselect['starmed'].values[i]],\
                     [np.abs(discselect[medval]-np.abs(discselect[width1])).values[i],np.abs(discselect[medval]+np.abs(discselect[width2])).values[i]],color='blue')

    if plotbar:
        barselect = mod.loc[mod['comp']=='bar']
        plt.scatter(barselect['starmed'],np.abs(barselect[medval]),marker=modsymbols[imod],edgecolor='black',facecolor='black',s=30.)
        for i in range(0,barselect['starmed'].values.size):
            plt.plot([barselect['starmed'].values[i],barselect['starmed'].values[i]],\
                     [np.abs(barselect[medval]-np.abs(barselect[width1])).values[i],np.abs(barselect[medval]+np.abs(barselect[width2])).values[i]],color='black')


    if plotknot:
        try:
            knotselect = mod.loc[mod['comp']=='knot']
            plt.scatter(knotselect['starmed'],np.abs(knotselect[medval]),marker=modsymbols[imod],edgecolor='red',facecolor='black',s=30.)
            for i in range(0,barselect['starmed'].values.size):
                plt.plot([knotselect['starmed'].values[i],knotselect['starmed'].values[i]],\
                         [np.abs(knotselect[medval]-np.abs(knotselect[width1])).values[i],np.abs(knotselect[medval]+np.abs(knotselect[width2])).values[i]],color='red')
        except:
            print('No knot!',imod)

#plt.axis([0.,5.,0.,1000.])

plt.xlabel('radius (kpc)')
plt.ylabel('$L_z$ angular momentum (kpc km/s)')
plt.ylabel('angle (degrees)')

plt.tight_layout()
plt.savefig('figures/distanceerror_diagnostic.png',dpi=300)
