import glob,os

import numpy as np

from .localtools import *

def read_mock_file(infile,nsamples = 1000,distance_uncertainty_factor = 1.,xsun=[-8.3,0.0,0.03],vsun=[11.1,244.24,7.25]):

    A = np.genfromtxt(infile,skip_header=1)

    Mock = dict()
    Mock['l']          = A[:,0]
    Mock['b']          = A[:,1]
    Mock['dist']       = A[:,2]
    Mock['edist']      = A[:,3]*distance_uncertainty_factor
    Mock['pml']        = A[:,4]
    Mock['epml']       = A[:,5]
    Mock['pmb']        = A[:,6]
    Mock['epmb']       = A[:,7]
    Mock['pmlpmbcorr'] = A[:,8]
    Mock['vlos']       = A[:,9]
    Mock['evlos']      = A[:,10]
    Mock['apogee']     = A[:,11]
    Mock['bulge']      = A[:,12]
    Mock['feh']        = A[:,13]
    Mock['alpha']      = A[:,14]
    Mock['age']        = A[:,15]
    Mock['index']      = A[:,16] # throw these as integers to be helpful...

    try:
        A = np.genfromtxt(infile,skip_header=1,dtype='S20')
        #print(A[:,17])
        Mock['apogee_id'] = A[:,17]

    except:
        Mock['apogee_id'] = Mock['index'].astype('S20')

    # make velocities
    theta  = (np.pi/2.)-Mock['b']*np.pi/180.
    phi    = Mock['l']*np.pi/180.


    vr     = Mock['vlos']
    vphi   = Mock['pml']*Mock['dist']*4.74
    vtheta = -Mock['pmb']*Mock['dist']*4.74
    v0     = spherical_to_cartesian_velocities(phi,theta,vr,vphi,vtheta)
    x0     = spherical_to_cartesian_positions(phi,theta,Mock['dist'])

    Mock['u'],Mock['v'],Mock['w'] = v0
    Mock['u'] += vsun[0]
    Mock['v'] += vsun[1]
    Mock['w'] += vsun[2]

    Mock['x'],Mock['y'],Mock['z'] = x0

    Mock['x'] += xsun[0]
    Mock['y'] += xsun[1]
    Mock['z'] += xsun[2]

    Mock['r'] = np.sqrt(Mock['x']**2 + Mock['y']**2 + Mock['z']**2)
    Mock['R'] = np.sqrt(Mock['x']**2 + Mock['y']**2)

    Mock['Lx'] = Mock['y']*Mock['w']-Mock['z']*Mock['v']
    Mock['Ly'] = Mock['z']*Mock['u']-Mock['x']*Mock['w']
    Mock['Lz'] = Mock['x']*Mock['v']-Mock['y']*Mock['u']

    Mock['L'] = np.sqrt(Mock['Lx']*Mock['Lx']+Mock['Ly']*Mock['Ly']+Mock['Lz']*Mock['Lz'])

    # Monte Carlo errors on the angular momentum (and also galactic radius)
    Mock['eR'] = np.zeros([Mock['Lx'].size,nsamples])
    Mock['eLx'] = np.zeros([Mock['Lx'].size,nsamples])
    Mock['eLy'] = np.zeros([Mock['Lx'].size,nsamples])
    Mock['eLz'] = np.zeros([Mock['Lx'].size,nsamples])

    for s in range(0,Mock['L'].size):
        vr     = Mock['vlos'][s] + np.random.normal(0.,Mock['evlos'][s],nsamples)

        derr   = Mock['dist'][s] + np.random.normal(0.,Mock['edist'][s],nsamples)
        vphi   = (Mock['pml'][s] + np.random.normal(0.,Mock['epml'][s],nsamples))*derr*4.74
        vtheta = -(Mock['pmb'][s] + np.random.normal(0.,Mock['epmb'][s],nsamples))*derr*4.74
        v0     = spherical_to_cartesian_velocities(phi[s],theta[s],vr,vphi,vtheta)
        x0     = spherical_to_cartesian_positions(phi[s],theta[s],derr)

        tmpu,tmpv,tmpw = v0
        tmpu += vsun[0]
        tmpv += vsun[1]
        tmpw += vsun[2]

        tmpx,tmpy,tmpz = x0
        tmpx += xsun[0]
        tmpy += xsun[1]
        tmpz += xsun[2]

        Mock['eR'][s]  = np.sqrt(tmpx*tmpx + tmpy*tmpy)
        Mock['eLx'][s] = tmpy*tmpw-tmpz*tmpv
        Mock['eLy'][s] = tmpz*tmpu-tmpx*tmpw
        Mock['eLz'][s] = tmpx*tmpv-tmpy*tmpu

    return Mock



def make_rotated_probabilities(Mock,good,COMPS,nanprint=False,nonorm=False):
    """
    Mock
    good
    COMPS: dictionary of component centroids

    """
    nparams = len(COMPS[0])

    probcomp = np.zeros([len(good),len(COMPS.keys())])

    for compnum,comp in enumerate(COMPS.keys()):

        # help the rotation go more smoothly
        ca = np.cos(COMPS[comp][7])
        sa = np.sin(COMPS[comp][7])

        # rotate the model components
        dlxm = COMPS[comp][1] * ca - COMPS[comp][2] * sa
        dlym = COMPS[comp][1] * sa + COMPS[comp][2] * ca
        dlzm = COMPS[comp][3]

        # rotate the hyperparameter model
        c11m = COMPS[comp][4]*ca**2. + COMPS[comp][5]*sa**2.
        c12m = COMPS[comp][4]*sa*ca  - COMPS[comp][5]*sa*ca
        c13m = 0.

        c21m = COMPS[comp][4]*sa*ca  - COMPS[comp][5]*sa*ca
        c22m = COMPS[comp][4]*sa**2. + COMPS[comp][5]*ca**2.
        c23m = 0.

        c31m = 0.
        c32m = 0.
        c33m = COMPS[comp][6]

        for inum,indx in enumerate(good):

            inputx = Mock['eLx'][indx]
            inputy = Mock['eLy'][indx]

            Cf = np.cov(np.array([inputx,inputy,Mock['eLz'][indx]])) # is this all I need?

            # add the hyperparameters
            Cf[0][0] += c11m
            Cf[0][1] += c12m
            Cf[1][0] += c21m
            Cf[1][1] += c22m
            Cf[2][2] += c33m

            # we are gettinng bad values here
            # one option is log determinant, e.g. https://numpy.org/doc/stable/reference/generated/numpy.linalg.slogdet.html#numpy.linalg.slogdet
            (signCf, logdetCf) = np.linalg.slogdet(Cf)
            covdet = np.abs(signCf * np.exp(logdetCf))
            covinv = np.linalg.inv(Cf)

            ldiff = np.array([Mock['Lx'][indx] - dlxm,
                              Mock['Ly'][indx] - dlym,
                              Mock['Lz'][indx] - dlzm])


            prefac = 1.
            # compute the un-normalised probability (eq 2)
            N1 = 1./np.sqrt(np.power(2.*np.pi,3.)*covdet)

            # put in a hard block check for very small values
            N2b = np.dot(ldiff.T,np.dot(covinv,ldiff))
            if N2b < 0:
                N2b = 20.
            N2 = np.exp(-prefac*N2b)

            # is this the minimum floor we should apply?
            if N1<1.e-12: N1 = 1.e-12
            if N2<1.e-12: N2 = 1.e-12

            N = N1*N2

            # check for bad component fits
            if (np.isnan(N)) & (nanprint==True):
                print("Compnum=",compnum,N,np.linalg.det(Cf),ldiff,np.exp(-np.dot(ldiff.T,np.dot(np.linalg.inv(Cf),ldiff))))
                print("Cf=",Cf)
                print("C=",C)

            probcomp[inum,compnum] = N

    # normalise for fitted probabilities:
    for compnum,comp in enumerate(COMPS.keys()):
        probcomp[:,compnum] *= COMPS[comp][0]

    # complete the normalisation on a per-star basis
    compsums = np.nansum(probcomp,axis=1)

    if not nonorm:
        for compnum,comp in enumerate(COMPS.keys()):
            probcomp[:,compnum] /= compsums

    # return array of probabilities
    return probcomp




def make_posterior_list_three_rotation_sorted(inputfile):
    """
    Reads a three-component posterior from the input file and computes component statistics.

    Args:
        inputfile (str): Path to the input file.

    Returns:
        tuple: A tuple containing component statistics (COMPS) and complete statistics (CStats).

    """
    A = np.genfromtxt(inputfile)


    CStats = dict()
    CStats['Likelihood'] = A[:,21]
    cats = ['tmpL','phi','theta','sxinv','syinv','szinv','alpha','f']



    for cnum in [0,1,2]:
        CStats[cnum] = dict()
        for icat,cat in enumerate(cats):
            CStats[cnum][cat] = A[:,cnum*len(cats)+icat]


    totalf = CStats[0]['f']+CStats[1]['f']+CStats[2]['f']
    CStats[0]['f'] /= totalf
    CStats[1]['f'] /= totalf
    CStats[2]['f'] /= totalf

    # make component products
    COMPS = dict()

    for cnum in range(0,3):

        if cnum==0:
            CStats[0]['logL'] = CStats[0]['tmpL']
            r   = 10.**CStats[0]['logL']

        if cnum==1:
            CStats[1]['logL'] = np.log10(CStats[1]['tmpL']*(10.**CStats[0]['logL']))
            r   = 10.**CStats[1]['logL']

        if cnum==2:
            CStats[2]['logL'] = np.log10(CStats[2]['tmpL']*CStats[1]['tmpL']*(10.**CStats[0]['logL']))
            r   = 10.**CStats[2]['logL']

        th  = np.arccos(CStats[cnum]['theta'])
        phi = CStats[cnum]['phi']
        x   = r*np.sin(th)*np.cos(phi)
        y   = r*np.sin(th)*np.sin(phi)
        z   = r*np.cos(th)
        CStats[cnum]['Lx'] = x
        CStats[cnum]['Ly'] = y
        CStats[cnum]['Lz'] = z

        COMPS[cnum] = [np.nanmedian(CStats[cnum]['f']),np.nanmedian(x),np.nanmedian(y),np.nanmedian(z),np.nanmedian(1./CStats[cnum]['sxinv']),np.nanmedian(1./CStats[cnum]['syinv']),np.nanmedian(1./CStats[cnum]['szinv']),np.nanmedian(CStats[cnum]['alpha'])]

    return COMPS,CStats




def cstats_to_comps(CStats, indx):
    """
    Converts cstats data to components data.

    Args:
        CStats (list): A list containing cstats data.
        indx (int): The index to extract data from CStats.

    Returns:
        dict: A dictionary containing components data.

    """
    COMPS = dict()

    for cnum in range(0, 3):
        r = 10. ** CStats[cnum]['logL']
        th = np.arccos(CStats[cnum]['theta'])
        phi = CStats[cnum]['phi']
        x = r * np.sin(th) * np.cos(phi)
        y = r * np.sin(th) * np.sin(phi)
        z = r * np.cos(th)

        if 'alpha' in CStats[cnum].keys():
            COMPS[cnum] = [
                CStats[cnum]['f'][indx],
                x[indx],
                y[indx],
                z[indx],
                1. / CStats[cnum]['sxinv'][indx],
                1. / CStats[cnum]['syinv'][indx],
                1. / CStats[cnum]['szinv'][indx],
                CStats[cnum]['alpha'][indx]
            ]
        else:
            COMPS[cnum] = [
                CStats[cnum]['f'][indx],
                x[indx],
                y[indx],
                z[indx],
                1. / CStats[cnum]['sxinv'][indx],
                1. / CStats[cnum]['syinv'][indx],
                1. / CStats[cnum]['szinv'][indx]
            ]

    return COMPS

def make_all_probabilities(data, criteria, CStats, nchains=100, nanprint=False, nonorm=False):
    """
    Computes all probabilities based on given data, criteria, and cstats.

    Args:
        data (dict): A dictionary containing data.
        criteria (ndarray): An array of criteria.
        CStats (list): A list containing cstats data.
        nchains (int, optional): The number of chains. Defaults to 100.

    Returns:
        tuple: A tuple containing all probabilities, percentile probabilities, and error probabilities.

    """
    # nchains = 100#CStats[0]['f'].size
    nstars = criteria.size
    print("Number of stars: {0}, with median R={1:4.3f} <sigma_Lz>={2:5.1f}".format(nstars, np.nanmedian(data['R'][criteria]),np.nanmedian(np.nanstd(data['eLz'][criteria],axis=1))))

    allprobs = np.zeros([nchains, nstars, 3])

    for indx in range(0, nchains):
        Cout = cstats_to_comps(CStats, indx)
        allprobs[indx] = make_rotated_probabilities(data, criteria, Cout,nanprint=nanprint,nonorm=nonorm)  # make_probabilities(AllDiscSNR, criteria, Cout)

    percentileprob = np.nanpercentile(allprobs, 50, axis=0)
    errorprob1 = np.nanpercentile(allprobs, 86, axis=0) - percentileprob
    errorprob2 = percentileprob - np.nanpercentile(allprobs, 14, axis=0)
    errorprob = np.nanmax([errorprob1,errorprob2],axis=0) # choose whichever side is more uncertain

    # probcomp = make_rotated_probabilities(AllDiscSNR,criteria,COMPS)
    return allprobs, percentileprob, errorprob




def print_classification(DModel,criteria,radii,percentileprob,errorprob,disccomp,barcomp,knotcomp,printdir,mockanalysis):

    minrad,maxrad = radii[0],radii[1]

    DModel['discprob'] = np.zeros(DModel['x'].size)
    DModel['ediscprob'] = np.zeros(DModel['x'].size)

    if disccomp>=0:
        DModel['discprob'][criteria] = percentileprob[:,disccomp]
        DModel['ediscprob'][criteria] = errorprob[:,disccomp]


    DModel['barprob'] = np.zeros(DModel['x'].size)
    DModel['ebarprob'] = np.zeros(DModel['x'].size)

    if barcomp>=0:
        DModel['barprob'][criteria] = percentileprob[:,barcomp]
        DModel['ebarprob'][criteria] = errorprob[:,barcomp]


    DModel['knotprob'] = np.zeros(DModel['x'].size)
    DModel['eknotprob'] = np.zeros(DModel['x'].size)

    if knotcomp>=0:
        DModel['knotprob'][criteria] = percentileprob[:,knotcomp]
        DModel['eknotprob'][criteria] = errorprob[:,knotcomp]


    # combine knot,bar
    #DModel['barprob'] = np.zeros(DModel['x'].size)
    #DModel['barprob'][criteria] = percentileprob[:,1]+percentileprob[:,0]

    #DModel['knotprob'] = np.zeros(DModel['x'].size)
    #DModel['knotprob'][criteria] = 0.

    f = open(printdir+'/3Component_AllFeHCutMembership_Percentiles_reduceSNR_r{}R{}_cyl.csv'.format(minrad,maxrad),'w')

    print('APOGEE_ID, P_knot, s_knot, P_bar, s_bar, P_disc, s_disc, X, Y, Z, Lx, Ly, Lz',file=f)

    for i in range(0,len(DModel['barprob'])):

        # print all regardless
        if (DModel['barprob'][i] + DModel['discprob'][i] + DModel['knotprob'][i]) > 0.7:
            if mockanalysis:
                if DModel['bulge'][i]==1.0:
                    pname = -1*DModel['index'][i]
                else:
                    pname = DModel['index'][i]
            else:
                pname = DModel['apogee_id'][i]
            print(pname,',',DModel['knotprob'][i],',',DModel['eknotprob'][i],',',DModel['barprob'][i],',',DModel['ebarprob'][i],',',DModel['discprob'][i],',',DModel['ediscprob'][i],',',\
                  DModel['x'][i],',',DModel['y'][i],',',DModel['z'][i],',',\
                        DModel['Lx'][i],',',DModel['Ly'][i],',',DModel['Lz'][i],file=f)

    f.close()




def print_h5_classification(DModel,criteria,radii,allprobs,disccomp,barcomp,knotcomp,printdir,mockanalysis):

    minrad,maxrad = radii[0],radii[1]

    DModel['discprob'] = np.zeros(DModel['x'].size)
    DModel['ediscprob'] = np.zeros(DModel['x'].size)

    if disccomp>=0:
        DModel['discprob'][criteria] = percentileprob[:,disccomp]
        DModel['ediscprob'][criteria] = errorprob[:,disccomp]


    DModel['barprob'] = np.zeros(DModel['x'].size)
    DModel['ebarprob'] = np.zeros(DModel['x'].size)

    if barcomp>=0:
        DModel['barprob'][criteria] = percentileprob[:,barcomp]
        DModel['ebarprob'][criteria] = errorprob[:,barcomp]


    DModel['knotprob'] = np.zeros(DModel['x'].size)
    DModel['eknotprob'] = np.zeros(DModel['x'].size)

    if knotcomp>=0:
        DModel['knotprob'][criteria] = percentileprob[:,knotcomp]
        DModel['eknotprob'][criteria] = errorprob[:,knotcomp]


    # combine knot,bar
    #DModel['barprob'] = np.zeros(DModel['x'].size)
    #DModel['barprob'][criteria] = percentileprob[:,1]+percentileprob[:,0]

    #DModel['knotprob'] = np.zeros(DModel['x'].size)
    #DModel['knotprob'][criteria] = 0.

    f = open(printdir+'/3Component_AllFeHCutMembership_Percentiles_reduceSNR_r{}R{}_cyl.csv'.format(minrad,maxrad),'w')

    print('APOGEE_ID, P_knot, s_knot, P_bar, s_bar, P_disc, s_disc, X, Y, Z, Lx, Ly, Lz',file=f)

    for i in range(0,len(DModel['barprob'])):
        #if (DModel['barprob'][i] > 0.7) | (DModel['discprob'][i] > 0.7) | (DModel['knotprob'][i] > 0.7):

        # print all regardless
        if (DModel['barprob'][i] + DModel['discprob'][i] + DModel['knotprob'][i]) > 0.7:
            if mockanalysis:
                if DModel['bulge'][i]==1.0:
                    pname = -1*DModel['index'][i]
                else:
                    pname = DModel['index'][i]
            else:
                pname = DModel['apogee_id'][i]
            print(pname,',',DModel['knotprob'][i],',',DModel['eknotprob'][i],',',DModel['barprob'][i],',',DModel['ebarprob'][i],',',DModel['discprob'][i],',',DModel['ediscprob'][i],',',\
                  DModel['x'][i],',',DModel['y'][i],',',DModel['z'][i],',',\
                        DModel['Lx'][i],',',DModel['Ly'][i],',',DModel['Lz'][i],file=f)

    f.close()
