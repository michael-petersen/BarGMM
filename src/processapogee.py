

# basic imports
import numpy as np
from numpy.linalg import eig, inv

# plotting elements
import matplotlib.pyplot as plt
%matplotlib inline
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

# to read fits, we need astropy
from astropy.io import fits

# transformation helpers
from localtools import *

# define the solar position relative to the galactic centre
usol,vsol,wsol = 9.,244.,9.
xsol,ysol,zsol = -8.178,0.,0.021

usol,vsol,wsol = 9.,244.,9.
xsol,ysol,zsol = -8.3,0.,0.021


filename = "final-apo17Gaiadr3-inner5kpc-sample.fits"

hdu = fits.open(filename)[1]
data = hdu.data
