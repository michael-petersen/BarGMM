
# basic imports
import numpy as np

# to read fits, we use astropy
from astropy.io import fits

# transformation helpers
from localtools import *

# define the solar position relative to the galactic centre (Jorge selected)
usol,vsol,wsol = 9.,244.,9.
xsol,ysol,zsol = -8.3,0.,0.021

indir = "data/"
basename = "APOGEE_all_feh8_SNRc"

filename = indir+"apogee/raw/APOGEE_BailerJones.fits"
hdu = fits.open(filename)[1]
data = hdu.data

# define a data quality mask:
# key parameters here are Fe/H cut and two quality cuts (defined correlations and ruwe)
good = np.isfinite(data['pmra']) &\
       np.isfinite(data['pmdec']) &\
       (data['FE_H']>-0.8) &\
       (data['pmra_pmdec_corr']<1.0) &\
       (data['ruwe']<1.4)


# create a data dictionary
APOGEE              = dict()
APOGEE['ra']        = data['ra'][good]
APOGEE['dec']       = data['dec'][good]
APOGEE['vlos']      = data['VHELIO_AVG'][good]
APOGEE['vscat']     = data['VSCATTER'][good]
APOGEE['verr']      = data['VERR'][good]
APOGEE['evlos']     = np.sqrt(data['VERR'][good]**2 + data['VSCATTER'][good])

# if using AstroNN distances
#APOGEE['dist']     = data['Dist'][good]*0.001 # let's do kpc
#APOGEE['edist']    = data['Dist_error'][good]*0.001 # let's do kpc

# if using Bailer-Jones distances
APOGEE['dist']      = data['GAIAEDR3_R_MED_PHOTOGEO'][good]*0.001 # let's do kpc
APOGEE['edist']     = 0.5*(data['GAIAEDR3_R_HI_PHOTOGEO'][good]-data['GAIAEDR3_R_LO_PHOTOGEO'][good])*0.001 # let's do kpc

APOGEE['pmra']      = data['pmra'][good]
APOGEE['epmra']     = data['pmra_error'][good]
APOGEE['pmdec']     = data['pmdec'][good]
APOGEE['epmdec']    = data['pmdec_error'][good]
APOGEE['pmcorr']    = data['pmra_pmdec_corr'][good]
APOGEE['ruwe']      = data['ruwe'][good]
APOGEE['index']     = np.arange(0,len(data['ruwe']),1)[good]
APOGEE['apogee_id'] = data['APOGEE_ID'][good]
APOGEE['feh']       = data['FE_H'][good]
APOGEE['alpha']     = data['MG_FE'][good]
APOGEE['age']       = np.zeros(APOGEE['ra'].size) - 99.0#data['age'][good]
APOGEE['apogee']    = np.ones(APOGEE['ra'].size)
APOGEE['comp']      = np.zeros(APOGEE['ra'].size) - 99.0

print('Accepted {} stars.'.format(APOGEE['comp'].size))


# convert to galactic coordinates
APOGEE['l'],APOGEE['b'] = rotate_galactic(APOGEE['ra']*np.pi/180.,APOGEE['dec']*np.pi/180.)
APOGEE['ldeg'],APOGEE['bdeg'] = APOGEE['l']*180./np.pi,APOGEE['b']*180./np.pi

# add cartesian to the dictionary
APOGEE['x'],APOGEE['y'],APOGEE['z'] = rotate_positions(APOGEE['ra']*np.pi/180.,APOGEE['dec']*np.pi/180.,APOGEE['dist'])
APOGEE['d'] = np.sqrt(APOGEE['x']*APOGEE['x']+APOGEE['y']*APOGEE['y']+APOGEE['z']*APOGEE['z'])

# apply solar offset
APOGEE['x'] += xsol
APOGEE['y'] += ysol
APOGEE['z'] += zsol
APOGEE['r'] = np.sqrt(APOGEE['x']*APOGEE['x']+APOGEE['y']*APOGEE['y']+APOGEE['z']*APOGEE['z'])
APOGEE['R'] = np.sqrt(APOGEE['x']*APOGEE['x']+APOGEE['y']*APOGEE['y'])

# define the stars we are going to use
criteria = APOGEE['R']<5.


# compute the velocities in galactic coordinates
APOGEE['pml'] = np.zeros(len(APOGEE['ra']))
APOGEE['pmb'] = np.zeros(len(APOGEE['ra']))
APOGEE['epml'] = np.zeros(len(APOGEE['ra']))
APOGEE['epmb'] = np.zeros(len(APOGEE['ra']))
APOGEE['pmlpmbcorr'] = np.zeros(len(APOGEE['ra']))

for i in range(0,len(APOGEE['ra'])):
    APOGEE['pml'][i],APOGEE['pmb'][i] = rotate_velocities(APOGEE['ra'][i]*np.pi/180.,APOGEE['dec'][i]*np.pi/180.,APOGEE['pmra'][i],APOGEE['pmdec'][i])
    cov = rotate_errors(APOGEE['ra'][i]*np.pi/180.,APOGEE['dec'][i]*np.pi/180.,APOGEE['epmra'][i],APOGEE['epmdec'][i],APOGEE['pmcorr'][i])

    # check for bad data
    if (cov[0][0] < 0.0) | (cov[1][1] < 0.0):
        print(i,APOGEE['epmra'][i],APOGEE['epmdec'][i],APOGEE['pmcorr'][i])

    APOGEE['epml'][i] = np.sqrt(cov[0][0])
    APOGEE['epmb'][i] = np.sqrt(cov[1][1])
    APOGEE['pmlpmbcorr'][i] = cov[0][1]/(APOGEE['epml'][i]*APOGEE['epmb'][i])



# make cartesian velocities
theta  = (np.pi/2.)-APOGEE['b']
phi    = APOGEE['l']
vr     = APOGEE['vlos']
vphi   = APOGEE['pml']*APOGEE['dist']*4.74
vtheta = -APOGEE['pmb']*APOGEE['dist']*4.74
v0     = spherical_to_cartesian_velocities(phi,theta,vr,vphi,vtheta)
x0     = spherical_to_cartesian_positions(phi,theta,APOGEE['dist'])

APOGEE['u'],APOGEE['v'],APOGEE['w'] = v0
APOGEE['u'] += usol
APOGEE['v'] += vsol
APOGEE['w'] += wsol

APOGEE['Lx'] = APOGEE['y']*APOGEE['w']-APOGEE['z']*APOGEE['v']
APOGEE['Ly'] = APOGEE['z']*APOGEE['u']-APOGEE['x']*APOGEE['w']
APOGEE['Lz'] = APOGEE['x']*APOGEE['v']-APOGEE['y']*APOGEE['u']

APOGEE['L'] = np.sqrt(APOGEE['Lx']*APOGEE['Lx']+APOGEE['Ly']*APOGEE['Ly']+APOGEE['Lz']*APOGEE['Lz'])

# this is the format for MULTINEST fitting
filename = indir+basename+".txt"
print_data(APOGEE,filename,criteria)

# this format includes the string names of the stars for cross-matching
filename = indir+basename+"_tagged.txt"
print_data_with_tags(APOGEE,filename,criteria)
