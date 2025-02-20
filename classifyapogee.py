
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
from models.apogee import *

print(modeltag+appendix)
Stars = read_mock_file(datafile)


classify = True


# specify which keys are being plotted
pltkeys = ['f','Lx', 'Ly','Lz','alpha','sxinv', 'syinv', 'szinv']


F = parameter_figure(pltkeys)

A = table_print(format='terminal')
B = table_print(format='markdown',outputfile=open(inputdir+'fits/README{0}.md'.format(modelappendix),'w'))
C = table_print(format='csv',outputfile=open(inputdir+'fits/table{0}.csv'.format(modelappendix),'w'))
D = table_print(format='tex',outputfile=open(inputdir+'fits/table{0}.tex'.format(modelappendix),'w'))

for p in [A,B,C,D]:
    p.print_header()

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
            F._add_data(comptag[minrad][cnum],ikey,radval,median,lo,hi)

        for p in [A,B,C,D]:
            p.print_columnX(comptag,minrad,cnum)

for p in [B,C,D]:
    p.f.close()


F._frame_figure()
F._print_figure(modeltag,appendix)

# generate classifications
if classify:
    # loop through all radii
    for irad,rads in enumerate(radii):

        # get the min/max cylindrical radii in kpc
        minrad,maxrad = rads[0],rads[1]

        # open the correct chain
        directory = inputdir+"fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
        inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
        COMPS,CStats = make_posterior_list_three_rotation_sorted(inputfile)

        # define the stars we want to classify
        #criteria = np.where((Stars['R']>(minrad)) & (Stars['R']<(maxrad)))[0]
        criteria = np.where((Stars['R']>(binprefacs[irad]*minrad)) & (Stars['R']<(binprefacs[irad]*maxrad)))[0]

        # how many random draws to take?
        nsamples = 500

        # do the probabilistic classification
        # allprobs is the full set of classifications
        # percentileprob is the 50h percentile
        # errorprob is the 1-sigma error on the classification
        allprobs,percentileprob,errorprob = make_all_probabilities(Stars,criteria,CStats,nchains=nsamples)


        # which component is which?
        disccomp,barcomp,knotcomp = compnum[minrad]

        if minrad==radii[0][0]:
            f = h5py.File(inputdir+"classifications/AllClassifications_{0}{1}_500.h5".format(modeltag,appendix),"w")

        for indx,starnum in enumerate(criteria):
            probabilities = np.zeros([nsamples,4]) # always disc, bar, knot, [x,y,z,Lx,Ly,Lz]
            if disccomp>=0: probabilities[:,0] = allprobs[:,indx,disccomp]
            if barcomp>=0:  probabilities[:,1] = allprobs[:,indx,barcomp]
            if knotcomp>=0: probabilities[:,2] = allprobs[:,indx,knotcomp]
            probabilities[0:12,3] = [Stars['R'][starnum],Stars['x'][starnum],Stars['y'][starnum],Stars['z'][starnum],Stars['Lx'][starnum],Stars['Ly'][starnum],Stars['Lz'][starnum],Stars['u'][starnum],Stars['v'][starnum],Stars['w'][starnum],np.nanstd(Stars['eR'][starnum]),np.nanstd(Stars['eLz'])]
            if binprefacs[irad] > 0.5:
                try:
                    dset = f.create_dataset(Stars['apogee_id'][starnum], data=probabilities)
                except:
                    print(Stars['apogee_id'][starnum])
            else:
                dset = f.create_dataset(Stars['apogee_id'][starnum]+b'*', data=probabilities)

        if minrad==radii[-1][0]: f.close()
