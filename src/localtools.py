
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


"""
Python PSP (Phase-Space Protocol) reader

MSP 25 Oct 2014 in original form
MSP  3 Dec 2015 committed to exptool
MSP  7 Mar 2016 constructed to theoretically handle niatr/ndatr
MSP 27 Aug 2016 added compatibility for dictionary support, the long-term goal of the reader once I commit to re-engineering everything.
MSP  8 Dec 2016 cleaned up subdividing inputs. needs much more cleaning, particularly eliminating many 'self' items from the Input class. Should also set up dictionary dump by default, could just engineer in at the end?
MSP 11 Mar 2019 set up to read yaml-derived input files. A method to diagnose problems would be amazing--currently written elsewhere.
MSP 14 Aug 2019 handle indexing=True from exp component inputs
MSP 17 Dec 2019 major revision to simplify
MSP 28 Sep 2021 deprecate parallelisms (move to particle.py)

PSP is a file format used by the EXP basis function expansion N-body code
written by Martin Weinberg.


TODO
-add protection for missing yaml?
-add handling for multiple components simultaneously
-add _make_dataframe

"""

import numpy as np

# requires yaml support: likely needs to be installed.
import yaml


class Input:
    """Python reader for Phase-Space Protocol (PSP) files used by EXP.

    inputs
    ----------
    filename : str
        The PSP filename to load.

    """

    def __init__(self, filename,comp=None, legacy=False,nbodies=-1,verbose=0,nout=-1):
        """
        inputs
        ------------
        comp    : str
            name of the component to return.
        legacy  : boolean
            if True, adds support for other exptool methods. unneeded if building from scratch.
        nbodies : integer
            reduce the number of bodies that are returned.
        verbose : integer
            levels of verbosity.
        nout    : integer
            deprecated compatibility parameter. use nbodies instead.

        """


        self.filename = filename
        self.nbodies = nbodies

        # deprecated.
        self.infile = filename

        # deprecated.
        self.nbodies = int(np.nanmax([nbodies,nout]))

        # initial check for file validity
        try:
            f = open(self.filename, 'rb')
            f.close()
        except Exception:
            raise IOError('Failed to load header from file "{}" - are you sure '
                          'this is a PSP file?'.format(filename))

        # test for split PSP files
        # TODO

        # do an initial read of the header
        self.header = self._read_primary_header()

        _comps = list(self.header.keys())


        # if a component is defined, retrieve data
        if comp != None:
            if comp not in _comps:
                raise IOError('The specified component does not exist.')

            else:
                self.data = self._read_component_data(self.header[comp])

        # if no header is defined, you will get just the primary header

        if (legacy) & (comp!=None):
            self.header = self.header[comp]
            self.comp = comp
            self._make_backward_compatible()
        elif legacy:
            raise IOError('A component must be specified for legacy usage.')


    def _read_component_header(self, f, comp_idx):
        """read in the header for a single component"""

        _ = f.tell()  # byte position of this component

        # TODO: if PSP changes, this will have to be altered
        if self._float_len == 4:
            _1,_2, nbodies, nint_attr, nfloat_attr, infostringlen = np.fromfile(f, dtype=np.uint32, count=6)

        else:
            nbodies, nint_attr, nfloat_attr, infostringlen = np.fromfile(f, dtype=np.uint32, count=4)

        # information string from the header
        head = np.fromfile(f, dtype=np.dtype((np.bytes_, infostringlen)),
                           count=1)

        # need the backward compatibility here.
        try:
            head_normal = head[0].decode()
            head_dict = yaml.safe_load(head_normal)
        except:
            # backward_compatibility
            head_dict = self._read_compatible_header(head)

        comp_data_pos = f.tell()  # byte position where component data begins

        # the default fields are (m, x, y, z, vx, vy, vz, p)
        nfields = 8
        comp_length = nbodies * (8 * int(head_dict['parameters']['indexing']) +
                                 self._float_len * nfields +
                                 4 * nint_attr +
                                 self._float_len * nfloat_attr)
        comp_data_end = f.tell() + comp_length  # byte pos. of comp. data end

        data = dict()
        data['index'] = comp_idx
        for k in head_dict:
            data[k] = head_dict[k]
        data['nint_attr'] = nint_attr
        data['nfloat_attr'] = nfloat_attr
        data['nbodies'] = nbodies
        data['data_start'] = comp_data_pos
        data['data_end'] = comp_data_end
        f.seek(comp_data_end)

        return data

    def _read_compatible_header(self,head):
        """read the old style of PSP header

        handling could be more general: this may have failure cases that I have not foreseen.

        """

        head_sep = head[0].decode().split(':')
        head_dict = dict()
        head_dict['parameters'] = dict()
        head_dict['parameters']['indexing'] = 0

        for istanza,stanza in enumerate(head_sep):

            if istanza==0:
                head_dict['name'] = stanza.strip()

            if istanza==1:
                head_dict['id'] = stanza.strip()

            if istanza > 1:
                stanza_sep = stanza.split(',')
                for param in stanza_sep:
                    head_dict['parameters'][param.split('=')[0].strip()] = param.split('=')[1].strip()

        return head_dict

    def _read_primary_header(self):
        """read the primary header of the PSP file"""

        primary_header = dict()
        nbodies = 0

        with open(self.filename, 'rb') as f:

            f.seek(16)  # find magic number
            cmagic, = np.fromfile(f, dtype=np.uint32, count=1)

            # check if it is float vs. double
            if cmagic == 2915019716:
                self._float_len = 4
                self._float_str = 'f'
            else:
                self._float_len = 8
                self._float_str = 'd'

            # reset to beginning and read current time
            f.seek(0)
            self.time, = np.fromfile(f, dtype='<f8', count=1)
            self._nbodies_tot, self._ncomp = np.fromfile(f, dtype=np.uint32,
                                                         count=2)

            for i in range(self._ncomp):
                data = self._read_component_header(f, i)
                primary_header[data.pop('name')] = data
                nbodies += data['nbodies']

        primary_header['nbodies'] = nbodies

        return primary_header

    def _read_component_data(self, comp_header):
        """read in all data for component"""

        dtype_str = []
        colnames = []
        if comp_header['parameters']['indexing']:
            # if indexing is on, the 0th column is Long
            dtype_str = dtype_str + ['l']
            colnames = colnames + ['index']

        dtype_str = dtype_str + [self._float_str] * 8
        colnames = colnames + ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'potE']

        dtype_str = dtype_str + ['i'] * comp_header['nint_attr']
        colnames = colnames + ['i_attr{}'.format(i)
                               for i in range(comp_header['nint_attr'])]

        dtype_str = dtype_str + [self._float_str] * comp_header['nfloat_attr']
        colnames = colnames + ['f_attr{}'.format(i)
                               for i in range(comp_header['nfloat_attr'])]

        dtype = np.dtype(','.join(dtype_str))

        out = np.memmap(self.filename,
                        dtype=dtype,
                        shape=(1, comp_header['nbodies']),
                        offset=int(comp_header['data_start']),
                        order='F', mode='r')

        tbl = dict()
        for i, name in enumerate(colnames):
            if self.nbodies > 0:
                tbl[name] = np.array(out['f{}'.format(i)][0], copy=True)[0:self.nbodies]
            else:
                tbl[name] = np.array(out['f{}'.format(i)][0], copy=True)

        del out  # close the memmap instance

        return tbl


    def _make_backward_compatible(self):
        """routine to make the dictionary style from above a drop-in replacement for old psp_io"""

        if self.nbodies > 0:
            self.mass = self.data['m'][0:self.nbodies]
            self.xpos = self.data['x'][0:self.nbodies]
            self.ypos = self.data['y'][0:self.nbodies]
            self.zpos = self.data['z'][0:self.nbodies]
            self.xvel = self.data['vx'][0:self.nbodies]
            self.yvel = self.data['vy'][0:self.nbodies]
            self.zvel = self.data['vz'][0:self.nbodies]
            self.pote = self.data['potE'][0:self.nbodies]

            try:
                self.indx = self.data['index'][0:self.nbodies]
            except:
                pass

        else:
            self.mass = self.data['m']
            self.xpos = self.data['x']
            self.ypos = self.data['y']
            self.zpos = self.data['z']
            self.xvel = self.data['vx']
            self.yvel = self.data['vy']
            self.zvel = self.data['vz']
            self.pote = self.data['potE']

            try:
                self.indx = self.data['index']
            except:
                pass

        # may also want to delete self.data in this case to save memory
        #del self.data

    def _make_dataframe(self):
        """routine to make the dictionary style from above pandas dataframe"""

        pass




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
    #print(cov)

    G = np.einsum('ab,ac->bc', pqgal,
                      np.einsum('ji,ik->jk', return_gaia_Agprime(), pqicrs))

    cov_to = np.einsum('ba,ac->bc', G,
                           np.einsum('ij,ki->jk', cov, G))

    return cov_to

#print(G)

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



def _make_backward_compatible(self):
    """routine to make the dictionary style from above a drop-in replacement for old psp_io"""

    self.mass = self.data['m']
    self.xpos = self.data['x']
    self.ypos = self.data['y']
    self.zpos = self.data['z']
    self.xvel = self.data['vx']
    self.yvel = self.data['vy']
    self.zvel = self.data['vz']
    self.pote = self.data['potE']

    try:
        self.indx = self.data['index']
    except:
        pass

    return self

    
