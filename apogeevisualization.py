

# basic imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# plotting utilities
import matplotlib.pyplot as plt;import matplotlib as mpl;import matplotlib.cm as cm;import matplotlib.colors as colors;from matplotlib import rc

majortickwidth,minortickwidth,fontsize = 1.5,0.75,10
majortickwidth,minortickwidth,fontsize = 1.0,0.5,10

cmap = mpl.cm.inferno # set a default perceptually-uniform colourmap
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['axes.linewidth'] = majortickwidth
for axname in ['xtick','ytick']:
    mpl.rcParams['{}.labelsize'.format(axname)] = fontsize
    mpl.rcParams['{}.major.width'.format(axname)] = majortickwidth
    mpl.rcParams['{}.minor.width'.format(axname)] = minortickwidth
    mpl.rcParams['{}.minor.visible'.format(axname)] = True



# data storage imports
import h5py

# local tools
from src.localtools import *
from src.fitmanagement import *

# select the model of interest
from models.apogee import *
F = pd.read_csv(inputdir+'fits/table.csv')
figname = ""


from models.bulgemock import *
F = pd.read_csv(inputdir+'fits/table_apog0.csv')
figname = "bulgemock"
"""
from models.barmock import *
F = pd.read_csv(inputdir+'fits/table_apog1.csv')
figname = "barmock"
"""

Stars = read_mock_file(datafile)


nsamples = 500

def gaussian(x, mu, sigma,A=1.0):
    return A*np.exp(-(x - mu)**2 / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))


bins = np.linspace(-500,1500,100)

xwidth = (bins[1]-bins[0])
xabscissa = bins[1:]-(0.5*xwidth)

fig = plt.figure(figsize=(7.5,10.5))
axdict = dict()

xmin,xmax,ymin,ymax = 0.09,0.95,0.05,0.99
nx,ny = 3,9
dx,dy = (xmax-xmin)/nx,(ymax-ymin)/ny

for j in range(0,ny):
    axdict[j] = dict()
    for i in range(0,nx):
        axdict[j][i] = fig.add_axes([xmin+i*dx,ymax-(j+1)*dy,dx,dy])


# loop through the components of the fits
#for index, row in F.iterrows():
for index, row in enumerate([x*0.5 for x in range(0,9)]):
    ttlgaussian = dict()

    lo,hi = row,row+1.0#row['binmin'],row['binmax']
    criteria = np.where((Stars['R']>lo)&(Stars['R']<hi))

    if (lo % 1.0) == 0:
        minrad,maxrad = int(lo),int(hi)
    else:
        minrad,maxrad = int(lo*10),int(hi*10)

    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    try:
        COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)
    except:
        continue

    # open a file to write all Gaussian parameters
    f = open(inputdir+'fits/allgaussians_{}.csv'.format(lo),'w')
    print('n,dim,mean1,width1,f1,mean2,width2,f2,mean3,width3,f3',file=f)

    medianeL = dict()
    medianeL['x'] = np.nanmedian(np.nanstd(Stars['eLx'][criteria],axis=1))
    medianeL['y'] = np.nanmedian(np.nanstd(Stars['eLy'][criteria],axis=1))
    medianeL['z'] = np.nanmedian(np.nanstd(Stars['eLz'][criteria],axis=1))

    for icomp,comp in enumerate(['x','y','z']):

        prefac = 1.0
        if comp=='z': prefac *= -1.0
        
        ttlgaussian[icomp] = np.zeros([nsamples,xabscissa.size])
        
        nstars = Stars['L{}'.format(comp)][criteria].size

        if comp=='z':
            hist, bin_edges = np.histogram(prefac*Stars['L{}'.format(comp)][criteria],bins=bins)
        elif comp=='x':
            hist, bin_edges = np.histogram(prefac*Stars['Ly'.format(comp)][criteria],bins=bins)
        elif comp=='y':
            hist, bin_edges = np.histogram(prefac*Stars['Lx'.format(comp)][criteria],bins=bins)

        inthist = np.nansum(hist)*xwidth
        axdict[index][icomp]

        for n in range(0,nsamples):
            print('{},{},'.format(n,comp),end="",file=f)
            ttlgaussian[icomp][n] = np.zeros(xabscissa.size)
            for key in [0,1,2]:
                #print(key,CStats[key][mean][n],np.sqrt(1./CStats[key]['{}inv'.format(wtag)][n]),CStats[key]['f'][n])
                
                gauss = gaussian(xabscissa,prefac*CStats[key]['L{}'.format(comp)][n],np.linalg.norm([np.sqrt(1./CStats[key]['s{}inv'.format(comp)][n]),medianeL[comp]]),A=CStats[key]['f'][n])
                
                if key!=2:
                    if comp=='x':
                        print('{},{},{},'.format(prefac*CStats[key]['Ly'.format(comp)][n],np.linalg.norm([np.sqrt(1./CStats[key]['syinv'.format(comp)][n]),medianeL[comp]]),CStats[key]['f'][n]),end="",file=f)
                    elif comp=='y':
                        print('{},{},{},'.format(prefac*CStats[key]['Lx'.format(comp)][n],np.linalg.norm([np.sqrt(1./CStats[key]['sxinv'.format(comp)][n]),medianeL[comp]]),CStats[key]['f'][n]),end="",file=f)
                    else:
                        print('{},{},{},'.format(prefac*CStats[key]['Lz'.format(comp)][n],np.linalg.norm([np.sqrt(1./CStats[key]['szinv'.format(comp)][n]),medianeL[comp]]),CStats[key]['f'][n]),end="",file=f)
                else:
                    if comp=='x':
                        print('{},{},{},'.format(prefac*CStats[key]['Ly'.format(comp)][n],np.linalg.norm([np.sqrt(1./CStats[key]['syinv'.format(comp)][n]),medianeL[comp]]),CStats[key]['f'][n]),file=f)
                    elif comp=='y':
                        print('{},{},{},'.format(prefac*CStats[key]['Lx'.format(comp)][n],np.linalg.norm([np.sqrt(1./CStats[key]['sxinv'.format(comp)][n]),medianeL[comp]]),CStats[key]['f'][n]),file=f)
                    else:
                        print('{},{},{},'.format(prefac*CStats[key]['Lz'.format(comp)][n],np.linalg.norm([np.sqrt(1./CStats[key]['szinv'.format(comp)][n]),medianeL[comp]]),CStats[key]['f'][n]),file=f)
                                # add a factor for potential rotation: 

                """
                if comp=='x': # x' = x cosa + y sina
                    gauss *= np.cos(CStats[key]['alpha'][n])
                    gauss2 = gaussian(xabscissa,CStats[key]['Ly'][n],np.linalg.norm([np.sqrt(1./CStats[key]['syinv'][n]),medianeL['y']]),A=CStats[key]['f'][n])
                    gauss += gauss2*np.sin(CStats[key]['alpha'][n])

                if comp=='y': # y' = y cosa - x sina
                    gauss *= np.cos(CStats[key]['alpha'][n])
                    gauss2 = gaussian(xabscissa,CStats[key]['Lx'][n],np.linalg.norm([np.sqrt(1./CStats[key]['sxinv'][n]),medianeL['x']]),A=CStats[key]['f'][n])
                    gauss += gauss2*np.sin(CStats[key]['alpha'][n])
"""
                
                axdict[index][icomp].plot(xabscissa,gauss,color='blue',lw=0.5,alpha=0.01)
                ttlgaussian[icomp][n] += gauss

        if (icomp==2):
            axdict[index][icomp].text(1.02,0.95,'{}$<$R$<${} kpc'.format(lo,hi),color='black',ha='left',va='top',rotation=90.,transform=axdict[index][key].transAxes)
        

    compname = ['Lx','Ly','Lz']
    for key in ttlgaussian.keys():
        for n in range(0,nsamples):
            axdict[index][key].plot(xabscissa,ttlgaussian[key][n],color='grey',lw=1.0,alpha=0.01)
        axdict[index][key].text(-400,0.0075,compname[key],color='black',ha='left',va='top')
        axdict[index][key].axis([np.nanmin(xabscissa),np.nanmax(xabscissa),0.,0.008])
        if (key!=0):
            axdict[index][key].set_yticklabels(())
        if ((key==0)&(index==0)):
            axdict[index][key].set_yticks([0.000,0.002,0.004,0.006,0.008])
        if ((key==0)&(index!=0)):
            axdict[index][key].set_yticks([0.000,0.002,0.004,0.006,0.008])
            axdict[index][key].set_yticklabels([0.000,0.002,0.004,0.006,''])
        ttlintegration = np.sum(ttlgaussian[key])
        #print(key,nstarsdict[key],ttlintegration)

axdict[8][1].set_xlabel('angular momentum (kpc km/s)')
axdict[0][0].set_ylabel('pdf')

plt.savefig('figures/fitcheck{}.png'.format(figname),dpi=300)
