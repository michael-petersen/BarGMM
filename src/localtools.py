
import numpy as np

def print_mock(catalog,filename,criteria1,used1,criteria2,used2):

    f = open(filename,'w')

    print("{0:10s} {1:10s} {2:7s} {3:7s} {4:7s} {5:7s} {6:7s} {7:7s} {8:7s} {9:9s} {10:9s} {11:7s} {12:7s} {13:7s} {14:7s} {15:7s} {16:7s}".format(
                       'l[deg]','b[deg]', 'dist[kpc]','edist[kpc]',
                       'pml[masyr]','epml[masyr]', 'pmb[masyr]', 'epmb[masyr]', 'pmlpmbcorr[]','vlos[kms]', 'evlos[kms]',
                       'apogee[0,1]', 'component[0,1]','feh[dex]', 'alpha[dex]', 'age[]', 'index'),file=f)

    nsources = len(catalog['l'][criteria1][used1])

    for i in range(0,nsources):
        print("{0:10.6f} {1:10.6f} {2:7.4f} {3:7.4f} {4:7.4f} {5:7.4f} {6:7.4f} {7:7.4f} {7:8.4f} {9:9.4f} {10:9.4f} {11:7.4f} {12:7.4f} {13:7.4f} {14:7.4f} {15:7.4f} {16:7.4f}".format(
                       catalog['l'][criteria1][used1][i],catalog['b'][criteria1][used1][i],
                       catalog['d'][criteria1][used1][i],catalog['edist'][criteria1][used1][i],
                       catalog['pml'][criteria1][used1][i],catalog['epml'][criteria1][used1][i],
                       catalog['pmb'][criteria1][used1][i],catalog['epmb'][criteria1][used1][i],
                       catalog['pmlpmbcorr'][criteria1][used1][i],
                       catalog['vlos'][criteria1][used1][i],catalog['evlos'][criteria1][used1][i],
                       catalog['apogee'][criteria1][used1][i],catalog['comp'][criteria1][used1][i],
                       catalog['feh'][criteria1][used1][i],catalog['alpha'][criteria1][used1][i],
                       catalog['age'][criteria1][used1][i],catalog['index'][criteria1][used1][i]),file=f)
        print("{0:10.6f} {1:10.6f} {2:7.4f} {3:7.4f} {4:7.4f} {5:7.4f} {6:7.4f} {7:7.4f} {7:8.4f} {9:9.4f} {10:9.4f} {11:7.4f} {12:7.4f} {13:7.4f} {14:7.4f} {15:7.4f} {16:7.4f}".format(
                       catalog['l'][criteria2][used2][i],catalog['b'][criteria2][used2][i],
                       catalog['d'][criteria2][used2][i],catalog['edist'][criteria2][used2][i],
                       catalog['pml'][criteria2][used2][i],catalog['epml'][criteria2][used2][i],
                       catalog['pmb'][criteria2][used2][i],catalog['epmb'][criteria2][used2][i],
                       catalog['pmlpmbcorr'][criteria2][used2][i],
                       catalog['vlos'][criteria2][used2][i],catalog['evlos'][criteria2][used2][i],
                       catalog['apogee'][criteria2][used2][i],catalog['comp'][criteria2][used2][i],
                       catalog['feh'][criteria2][used2][i],catalog['alpha'][criteria2][used2][i],
                       catalog['age'][criteria2][used2][i],catalog['index'][criteria2][used2][i]),file=f)

    f.close()





from scipy.interpolate import interp1d

def draw_sources(radk,numk,radj,nsources = 15000,TOL=0.03):
    """ draw sources from an input distribution to match a specified distribution

    inputs
    ---------
    radk     :
    numk     :
    radj     :
    nsources :

    """
    # interpolate the input distribution

    spl = interp1d(numk,radk,'nearest')

    #used = np.zeros(nsources,dtype='int')
    used = []
    desdist = np.zeros(nsources)

    for i in range(0,nsources):
        #if i%1500==0:print(i)
        r = np.random.rand()
        desired_dist = spl(r)
        desdist[i] = desired_dist

        # set up a distance tolerance of 1%.
        # should this be converted to a heliocentric distance?
        w = np.where( (radj > (1.0-TOL)*desired_dist) & (radj < (1.0+TOL)*desired_dist))[0]

        # randomly select a star from the list
        t = np.random.choice(w)

        draws = 0

        # allow duplicates?
        while t in used:
            t = np.random.choice(w)
            draws+=1
            if draws==w.size:
                #print('failure',w.size,desired_dist)
                t = -1#np.random.choice(w)
                break
                #t = 1000000+i

        #used[i] = t
        if t>0:
            used.append(t)

    print('Target num={}, returned num={}'.format(nsources,len(np.array(used))))
    return np.array(used),desdist,spl





def to_galactic(x0,y0,z0,u0,v0,w0,twopi=True):

    rad = np.sqrt(x0**2+y0**2+z0**2)
    xphi= np.arctan2(y0,x0)
    xth = np.arccos(z0/rad)

    xur = np.zeros([3,x0.size])
    xur[0]= np.sin(xth)*np.cos(xphi)
    xur[1]= np.sin(xth)*np.sin(xphi)
    xur[2]= np.cos(xth)

    xuth = np.zeros([3,x0.size])
    xuth[0]= np.cos(xth)*np.cos(xphi)
    xuth[1]= np.cos(xth)*np.sin(xphi)
    xuth[2]=-np.sin(xth)

    xuphi = np.zeros([3,x0.size])
    xuphi[0]=-np.sin(xphi)
    xuphi[1]=+np.cos(xphi)
    xuphi[2]= 0.

    vr =    u0*  xur[0] + v0*  xur[1] + w0*  xur[2]
    vth=    u0* xuth[0] + v0* xuth[1] + w0* xuth[2]
    vphi=   u0*xuphi[0] + v0*xuphi[1] + w0*xuphi[2]

    vb= -vth

    # match the astropy output
    vl= vphi

    dk  =4.74057           #conversion from km/s
    par =1./rad             #arc sec --> rad in [kpc]
    dmul=vl / dk * par
    dmub=vb / dk * par

    f=np.pi/180.
    dB=np.arcsin(z0/rad)/f

    if twopi:
        dL=np.arctan(y0/x0)/f
        #dL[(y0<0)&(x0>0.)] += 360.
        #dL[(y0>0)&(x0<0.)] += 180.
        #dL[(y0<0)&(x0<0.)] += 180.
    else:
        dL = np.arctan2(y0,x0)/f


    return dL,dB,rad,vr,dmul,dmub



# constants we will need
rad_to_deg = 180./np.pi      # radians to degrees
deg_to_rad = np.pi/180.      # degrees to radians
pc_to_km = 3.086e13          # parsec to km
s_to_myr = 3.17098e-14       # seconds to Myr
astronomicalG = 0.0043009125 # pc/Msun/km/s/km/s
mas_to_kms = 4.74            # scaling for proper motions to velocities (this can be more accurate??)

# do we need anything else?

# angular conversion definitions (helpful below)
def phi_theta(x,y,z):
  return (180/np.pi)*np.arctan2(y,x),(180/np.pi)*np.arctan(z/np.sqrt(x*x+y*y))

def spherical_to_cartesian_positions(phi,theta,r):
  cosp,sinp,cost,sint = np.cos(phi),np.sin(phi),np.cos(theta),np.sin(theta)
  x = r*sint*cosp
  y = r*sint*sinp
  z = r*cost
  return np.array([x,y,z])

def spherical_to_cartesian_velocities(phi,theta,vr,vphi,vtheta):
  cosp,sinp,cost,sint = np.cos(phi),np.sin(phi),np.cos(theta),np.sin(theta)
  vx = sint*cosp*vr + cost*cosp*vtheta - sinp*vphi
  vy = sint*sinp*vr + cost*sinp*vtheta + cosp*vphi
  vz = cost     *vr - sint     *vtheta
  return np.array([vx,vy,vz])




def return_gaia_Agprime():
    """return the matrix in eq 3.61, key to transform from ICRS to galactic coordinates"""
    return np.array([[-0.0548755604162154,-0.8734370902348850,-0.4838350155487132],
                     [+0.4941094278755837,-0.4448296299600112,+0.7469822444972189],
                     [-0.8676661490190047,-0.1980763734312015,+0.4559837761750669]])

def return_dummy():
    """return the matrix in eq 3.61, key to transform from ICRS to galactic coordinates"""
    return np.array([[1.0,0.0,0.0],
                     [0.0,1.0,0.0],
                     [0.0,0.0,1.0]])



def return_ricrs(a,d):
    """ eq. 3.57"""
    return np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)]).T

def return_picrs(a,d):
    """ eq. 3.64, unit vector of increasing alpha"""
    return np.array([-np.sin(a),np.cos(a),0.]).T

def return_qicrs(a,d):
    """ eq. 3.64, unit vector of increasing delta"""
    return np.array([-np.cos(a)*np.sin(d),-np.sin(a)*np.sin(d),np.cos(d)]).T

def return_muicrs(a,d,mua,mud):
    """ eq. 3.66, the proper motion vector"""
    p = return_picrs(a,d)
    q = return_qicrs(a,d)
    return np.dot(p,mua) + np.dot(q,mud)


def return_rgal(l,b):
    """ eq. 3.58"""
    return np.array([np.cos(l)*np.cos(b),np.sin(l)*np.cos(b),np.sin(b)]).T

def return_pgal(l,b):
    """ eq. 3.66, unit vector of increasing alpha"""
    return np.array([-np.sin(l),np.cos(l),0.]).T

def return_qgal(l,b):
    """ eq. 3.66, unit vector of increasing delta"""
    return np.array([-np.cos(l)*np.sin(b),-np.sin(l)*np.sin(b),np.cos(b)]).T

def return_mugal(l,b,mul,mub):
    """ eq. 3.66, the proper motion vector"""
    p = return_pgal(l,b)
    q = return_qgal(l,b)
    return np.dot(p,mul) + np.dot(q,mub)


def rotate_velocities(a,d,mua,mud,rotmat=return_gaia_Agprime()):
    """eq 3.68, """
    mu = return_muicrs(a,d,mua,mud)
    mugal = np.dot(rotmat,mu) # eq. 3.68

    # solve for positions
    ricrs = return_ricrs(a,d)
    rgal = np.dot(rotmat,ricrs)

    # implement eq 3.63
    ell,b = np.arctan2(rgal[1],rgal[0]),np.arctan2(rgal[2],np.sqrt(rgal[0]*rgal[0]+rgal[1]*rgal[1]))

    p = return_pgal(ell,b)
    q = return_qgal(ell,b)

    mul = np.dot(p.T,mugal)
    mub = np.dot(q.T,mugal)
    #print(mul,mub)
    return mul,mub



def rotate_errors(a,d,pmra_e,pmdec_e,pmcorr,rotmat=return_gaia_Agprime()):
    ricrs = return_ricrs(a,d)
    picrs = return_picrs(a,d)
    qicrs = return_qicrs(a,d)

    rgal = np.dot(rotmat,ricrs)

    # implement eq 3.63
    ell = np.arctan2(rgal[1],rgal[0])
    b = np.arctan2(rgal[2],np.sqrt(rgal[0]*rgal[0]+rgal[1]*rgal[1]))

    pgal = return_pgal(ell,b)
    qgal = return_qgal(ell,b)

    pqgal = np.stack((pgal, qgal), axis=-1)
    pqicrs = np.stack((picrs, qicrs), axis=-1)

    cov = np.array([[pmra_e*pmra_e,pmra_e*pmdec_e*pmcorr],[pmra_e*pmdec_e*pmcorr,pmdec_e*pmdec_e]])

    G = np.einsum('ab,ac->bc', pqgal,
                      np.einsum('ji,ik->jk', return_gaia_Agprime(), pqicrs))

    cov_to = np.einsum('ba,ac->bc', G,
                           np.einsum('ij,ki->jk', cov, G))

    return cov_to


def rotate_positions(a,d,dist,rotmat=return_gaia_Agprime()):
    """eq 3.68, """
    # solve for positions
    ricrs = return_ricrs(a,d).T
    rgal = np.dot(rotmat,ricrs)
    return dist*rgal


def rotate_galactic(a,d,rotmat=return_gaia_Agprime()):
    ricrs = return_ricrs(a,d).T

    rgal = np.dot(rotmat,ricrs)

    # implement eq 3.63
    ell = np.arctan2(rgal[1],rgal[0])
    b = np.arctan2(rgal[2],np.sqrt(rgal[0]*rgal[0]+rgal[1]*rgal[1]))

    return ell,b




def print_data(catalog,filename,criteria1):

    f = open(filename,'w')

    print("{0:10s} {1:10s} {2:7s} {3:7s} {4:7s} {5:7s} {6:7s} {7:7s} {8:7s} {9:9s} {10:9s} {11:7s} {12:7s} {13:7s} {14:7s} {15:7s} {16:7s}".format(
                       'l[deg]','b[deg]', 'dist[kpc]','edist[kpc]',
                       'pml[masyr]','epml[masyr]', 'pmb[masyr]', 'epmb[masyr]', 'pmlpmbcorr[]','vlos[kms]', 'evlos[kms]',
                       'apogee[0,1]', 'component[0,1]','feh[dex]', 'alpha[dex]', 'age[]', 'index'),file=f)

    nsources = len(catalog['l'][criteria1])

    for i in range(0,nsources):
        print("{0:10.6f} {1:10.6f} {2:7.4f} {3:7.4f} {4:7.4f} {5:7.4f} {6:7.4f} {7:7.4f} {8:7.4f} {9:9.4f} {10:9.4f} {11:7.4f} {12:7.4f} {13:7.4f} {14:7.4f} {15:7.4f} {16:7.4f}".format(
                       catalog['ldeg'][criteria1][i],catalog['bdeg'][criteria1][i],
                       catalog['d'][criteria1][i],catalog['edist'][criteria1][i],
                       catalog['pml'][criteria1][i],catalog['epml'][criteria1][i],
                       catalog['pmb'][criteria1][i],catalog['epmb'][criteria1][i],
                       catalog['pmlpmbcorr'][criteria1][i],
                       catalog['vlos'][criteria1][i],catalog['evlos'][criteria1][i],
                       catalog['apogee'][criteria1][i],catalog['comp'][criteria1][i],
                       catalog['feh'][criteria1][i],catalog['alpha'][criteria1][i],
                       catalog['age'][criteria1][i],catalog['index'][criteria1][i]),file=f)

    f.close()



def print_data_with_tags(catalog,filename,criteria1):

    f = open(filename,'w')

    print("{0:10s} {1:10s} {2:7s} {3:7s} {4:7s} {5:7s} {6:7s} {7:7s} {8:7s} {9:9s} {10:9s} {11:7s} {12:7s} {13:7s} {14:7s} {15:7s} {16:7s} {17:20s}".format(
                       'l[deg]','b[deg]', 'dist[kpc]','edist[kpc]',
                       'pml[masyr]','epml[masyr]', 'pmb[masyr]', 'epmb[masyr]', 'pmlpmbcorr[]','vlos[kms]', 'evlos[kms]',
                       'apogee[0,1]', 'component[0,1]','feh[dex]', 'alpha[dex]', 'age[]', 'index','apogeeID'),file=f)

    nsources = len(catalog['l'][criteria1])

    for i in range(0,nsources):
        print("{0:10.6f} {1:10.6f} {2:7.4f} {3:7.4f} {4:7.4f} {5:7.4f} {6:7.4f} {7:7.4f} {8:7.4f} {9:9.4f} {10:9.4f} {11:7.4f} {12:7.4f} {13:7.4f} {14:7.4f} {15:7.4f} {16:7.4f} {17:20s}".format(
                       catalog['ldeg'][criteria1][i],catalog['bdeg'][criteria1][i],
                       catalog['d'][criteria1][i],catalog['edist'][criteria1][i],
                       catalog['pml'][criteria1][i],catalog['epml'][criteria1][i],
                       catalog['pmb'][criteria1][i],catalog['epmb'][criteria1][i],
                       catalog['pmlpmbcorr'][criteria1][i],
                       catalog['vlos'][criteria1][i],catalog['evlos'][criteria1][i],
                       catalog['apogee'][criteria1][i],catalog['comp'][criteria1][i],
                       catalog['feh'][criteria1][i],catalog['alpha'][criteria1][i],
                       catalog['age'][criteria1][i],catalog['index'][criteria1][i],catalog['apogee_id'][criteria][i]),file=f)

    f.close()
