
import numpy as np

inputdir  = 'data/bulgemock/'
datafile  = inputdir+"ModelBarYoungBulgeMock10000.txt"
modeltag  = "ModelBarYoungBulgeMock10000"
modelname = 'BulgeMock'
mockanalysis = True

# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['disc','knot','knot']
comptag[1] = ['disc','knot','bar']
comptag[2] = ['disc','-','bar']
comptag[3] = ['-','-','disc']
comptag[4] = ['-','disc','-']

# by component number, which component is [disc,bar,knot]?
compnum = dict()
compnum[0] = [0,-1,2]
compnum[1] = [0,2,1]
compnum[2] = [0,2,-1]
compnum[3] = [2,-1,-1]
compnum[4] = [1,-1,-1]

# for ellipse tracing
# by [bar,disc,knot], which component is which?
complist = dict()
complist[0] = [1,2,2]
complist[1] = [1,2,0]
complist[2] = [1,-1,0]
complist[3] = [-1,-1,1]
complist[4] = [-1,1,-1]



from exptool.analysis import pattern
from exptool.analysis import trapping

indirtrap = '/Users/mpetersen/Notebooks/Dynamics/FixedPotential/data/Disk001/'

# polar trapping criteria
criteria = {}
criteria[0] = [(0.0,np.pi/6.),(33.,251.0),(0.0,0.001),(0.0,np.pi/16.),(0.,0.025)] # x1 orbits
criteria[1] = [(3*np.pi/8.,np.pi/2.),(33.,251.0),(0.0,0.001),(0.0,np.pi/8.),(0.,1.)] # x2 orbits
criteria[2] = [(0.0,np.pi/6.),(33.,251.0),(0.0,0.0010),(np.pi/16.,np.pi/8.),(0.,1.)] # generic bar supporting
#criteria[3] = [(0.0,np.pi/2.),(251.,10000.0),(0.0,1.0),(0.0,2.0),(0.,1.)] # nyquist excluded
criteria[4] = [(0.0,np.pi/6.),(33.,251.0),(0.001,0.0015),(np.pi/16.,np.pi/8.),(0.,1.)] # loose generic bar supporting
#criteria[5] = [(0.0,np.pi/2.),(10.,251.0),(0.00,0.005),(0,np.pi/16.),(0.02,1.)] # lagrange
#criteria[6] = [(0.0,np.pi/16.),(10.,251.0),(0.00,0.001),(np.pi/16.,np.pi/6.),(0.00,0.015)] # lagrange

trapping_array = {}
for critnum in criteria.keys():
    bar_times,trapping_array[critnum] = trapping.read_trapping_file(indirtrap+'trapping_criteria_080618_{}.dat'.format(critnum),tdtype='i1')
