
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
#Stars = read_mock_file(datafile)


#classify = True



xmin,dx = 0.1,0.16
ymin,dy = 0.1,0.16
xbuf,ybuf = 0.01,0.01

fig = plt.figure(figsize=(12,12),facecolor='white')

#for irad,rads in enumerate():
for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    print(minrad,maxrad)

    if minrad == 15:
        pass
    else:
        continue

    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
    A = np.genfromtxt(inputfile)
    CStats = dict()
    CStats['Likelihood'] = A[:,24]
    # the columns we want are 0,7,8,15,16,23
    nums = [0,7,8,15,16,23]
    tags = ['$\\log L$','$f_1$','$\\xi_2$','$f_2$','$\\xi_3$','$f_3$']

    for icat1,cat1 in enumerate(nums):
        for icat2,cat2 in enumerate(nums):
            if icat1<icat2:
                ax = fig.add_axes([xmin+icat1*(dx+xbuf),ymin+(icat2-1)*(dy+ybuf),dx,dy])
                for c in [0,1,2]:
                    xval = A[:,cat1]
                    yval = A[:,cat2]

                    ax.scatter(xval,yval,facecolor='black',edgecolor='None',s=1.,alpha=0.2)
                    if icat1==0:
                        ax.set_ylabel(tags[icat2])
                    if icat1+1==icat2:
                        ax.set_xlabel(tags[icat1])
                    if icat1+1<icat2:
                        ax.set_xticklabels(())
                    if icat1>0:
                        ax.set_yticklabels(())
                    ax.axis([0.,1.,0.,1.])

                    if icat1==0:
                        ax.axis([0.,3.,0.,1.])
                ax.tick_params(axis="both",direction="in",which="both")

plt.savefig('/Users/mpetersen/Downloads/testcorner.png',dpi=300)



"""

# specify which keys are being plotted
pltkeys = ['f','Lx', 'Ly','Lz','alpha','sxinv', 'syinv', 'szinv']


xmin,dx = 0.1,0.1
ymin,dy = 0.1,0.1
xbuf,ybuf = 0.01,0.01

fig = plt.figure(figsize=(12,12),facecolor='white')

#for irad,rads in enumerate():
for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    print(minrad,maxrad)

    if minrad == 25:
        pass
    else:
        continue

    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
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

    cvals = ['blue','red','black']
    for icat1,cat1 in enumerate(cats):
        for icat2,cat2 in enumerate(cats):
            if icat1<icat2:
                ax = fig.add_axes([xmin+icat1*(dx+xbuf),ymin+icat2*(dy+ybuf),dx,dy])
                for c in [0,1,2]:
                    xval = CStats[c][cat1]
                    if ('inv' in cat1):
                        xval = np.sqrt(1./CStats[c][cat1])
                    if (cat1 == 'tmpL'):
                        xval = np.log10(CStats[c][cat1])

                    yval = CStats[c][cat2]
                    if ('inv' in cat2):
                        yval = np.sqrt(1./CStats[c][cat2])
                    if (cat2 == 'tmpL'):
                        tval = np.log10(CStats[c][cat2])
                    ax.scatter(xval,yval,facecolor=cvals[c],edgecolor='None',s=1.)
                    if icat1==0:
                        ax.set_ylabel(cat2)
                    if icat1+1==icat2:
                        ax.set_xlabel(cat1)
                    if icat1+1<icat2:
                        ax.set_xticklabels(())
                    if icat1>0:
                        ax.set_yticklabels(())
                ax.tick_params(axis="both",direction="in",which="both")

plt.savefig('/Users/mpetersen/Downloads/testcorner.png',dpi=300)



fig = plt.figure(figsize=(15,14),facecolor='white')

xmin,dx = 0.05,0.25
ymin,dy = 0.05,0.25
xbuf,ybuf = 0.06,0.06
ax2 = fig.add_axes([xmin+3*(dx+xbuf)-xbuf+0.005,ymin+0*(dy+ybuf),0.01,dy])

#for irad,rads in enumerate():
i = 0
for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    print(minrad,maxrad)

    ax = fig.add_axes([xmin+int(np.floor(i/3))*(dx+xbuf),ymin+(2-int(np.floor(i%3)))*(dy+ybuf),dx,dy])
    i+=1
    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
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

    corrmat = np.zeros([24,24])
    for icat1 in range(0,24):
        for icat2 in range(icat1+1,24):
            corrmat[icat1,icat2] = np.corrcoef(np.array([A[:,icat1],A[:,icat2]]))[0][1]

    ax.imshow(np.abs(corrmat),vmin=0.,vmax=0.5,cmap=cm.Grays,extent=[0,24,0,24],aspect='auto')

    # add the components themselves
    disccomp,barcomp,knotcomp = compnum[minrad]

    cwheel = ['']
    # disc is limegreen
    # bar in firebrick
    # knot is indigo
    dcolor = 'limegreen'
    bcolor = 'firebrick'
    kcolor = 'indigo'

    # outline non-real components
    #ax.plot([16,16,24,24,16],[0,8,8,0,0],color='black',lw=3.)
    #ax.plot([8,8,16,16,8],[8,16,16,8,8],color='black',lw=3.)
    #ax.plot([0,0,8,8,0],[16,24,24,16,16],color='black',lw=3.)

    # overwrite squares with real components
    if disccomp==2:
        ax.plot([16,16,24,24,16],[0,8,8,0,0],color=dcolor,lw=3.)
    elif disccomp==1:
        ax.plot([8,8,16,16,8],[8,16,16,8,8],color=dcolor,lw=3.)
    elif disccomp==0:
        ax.plot([0,0,8,8,0],[16,24,24,16,16],color=dcolor,lw=3.)

    if barcomp==2:
        ax.plot([16,16,24,24,16],[0,8,8,0,0],color=bcolor,lw=3.)
    elif barcomp==1:
        ax.plot([8,8,16,16,8],[8,16,16,8,8],color=bcolor,lw=3.)
    elif barcomp==0:
        ax.plot([0,0,8,8,0],[16,24,24,16,16],color=bcolor,lw=3.)\

    if knotcomp==2:
        ax.plot([16,16,24,24,16],[0,8,8,0,0],color=kcolor,lw=3.)
    elif knotcomp==1:
        ax.plot([8,8,16,16,8],[8,16,16,8,8],color=kcolor,lw=3.)
    elif knotcomp==0:
        ax.plot([0,0,8,8,0],[16,24,24,16,16],color=kcolor,lw=3.)



    ax.set_xticks(np.arange(0,24,1)+0.5)
    ax.set_yticks(np.arange(0,24,1)+0.5)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_title('{}-{}kpc'.format(minrad*binprefacs[irad],maxrad*binprefacs[irad]))


    cats = ['$\\log |L|$','$\\phi_1$','$\\theta_1$','$\\sigma_{x,1}$','$\\sigma_{y,1}$','$\\sigma_{z,1}$','$\\alpha_1$','$f_1$',\
            '$\\xi_2$','$\\phi_2$','$\\theta_2$','$\\sigma_{x,2}$','$\\sigma_{y,2}$','$\\sigma_{z,2}$','$\\alpha_2$','$f_2$',\
            '$\\xi_3$','$\\phi_3$','$\\theta_3$','$\\sigma_{x,3}$','$\\sigma_{y,3}$','$\\sigma_{z,3}$','$\\alpha_3$','$f_3$']

    ax.set_yticklabels(cats[::-1])
    ax.set_xticklabels(cats,rotation=90)

cmin,cmax = 0.,0.5
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm.Grays,norm=norm)
cb1.set_label('correlation')
cb1.set_ticks([0.,0.25,0.5])
cb1.ax.minorticks_off()

plt.savefig('/Users/mpetersen/Downloads/testcorr.png',dpi=300)
"""


"""
fig = plt.figure(figsize=(5,4),facecolor='white')
ax = fig.add_axes([0.15,0.03,0.7,0.8])
ax2 = fig.add_axes([0.86,0.03,0.02,0.8])

#for irad,rads in enumerate():
for irad,rads in enumerate(radii):
    minrad,maxrad = rads[0],rads[1]
    print(minrad,maxrad)

    if minrad == 4:
        pass
    else:
        continue

    directory = inputdir+"/fits/{0}_d{1:02d}{2:02d}{3}".format(modeltag,minrad,maxrad,appendix)
    inputfile = directory+'/chains/gaussian-post_equal_weights.dat'
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

    cvals = ['blue','red','black']
    corrmat = np.zeros([24,24])
    for icat1 in range(0,24):
        for icat2 in range(icat1+1,24):
            corrmat[icat1,icat2] = np.corrcoef(np.array([A[:,icat1],A[:,icat2]]))[0][1]

ax.imshow(np.abs(corrmat),vmin=0.,vmax=0.5,cmap=cm.Grays,extent=[0,24,0,24])
ax.set_xticks(np.arange(0,24,1)+0.5)
ax.set_yticks(np.arange(0,24,1)+0.5)
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

cats = ['logL','phi0','theta0','sxinv0','syinv0','szinv0','alpha0','f0',\
        'fracL1','phi1','theta1','sxinv1','syinv1','szinv1','alpha1','f1',\
        'fracL2','phi2','theta2','sxinv2','syinv2','szinv2','alpha2','f2']

ax.set_yticklabels(cats[::-1])
ax.set_xticklabels(cats,rotation=90)

cmin,cmax = 0.,0.5
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm.Grays,norm=norm)
cb1.set_label('correlation')
cb1.set_ticks([0.,0.25,0.5])
cb1.ax.minorticks_off()

plt.savefig('/Users/mpetersen/Downloads/testcorr.png',dpi=300)
"""
