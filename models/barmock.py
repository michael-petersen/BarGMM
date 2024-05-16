
import numpy as np

inputdir     = 'data/barmock/'
datafile     = inputdir+"ModelBarYoungMock10000b.txt"
modeltag     = "ModelBarYoungMock10000b"
appendix     = "_apog0"  # if this is for the all-sky sample
appendix     = "_apog1"  # if this is for the SDSS-coverage sample

modelname    = 'BarMock'
mockanalysis = True

# what are the bins available to classify?
radii      = [[0,1],[1,2],[2,3],[3,4]]
binprefacs = [ 1.,   1.,   1.,   1.]


if appendix=="_apog0":
    apogeeflag = 0

    # what are the (hand-defined) components?
    comptag = dict()
    comptag[0] = ['disc','bar','knot']
    comptag[1] = ['disc','bar','knot']
    comptag[2] = ['disc','bar','knot']
    comptag[3] = ['disc','-','bar']

    # by component number, which component is [disc,bar,knot]?
    # assign the number of the cluster
    compnum = dict()
    compnum[0] = [0,1,2]
    compnum[1] = [0,1,2]
    compnum[2] = [0,1,2]
    compnum[3] = [0,2,-1]

    # for ellipse tracing
    # by [bar,disc,knot]=[0,1,2], which component is which?
    # assign the number of the component
    complist = dict()
    complist[0] = [1,0,2]
    complist[1] = [1,0,2]
    complist[2] = [1,0,2]
    complist[3] = [1,-1,0]

if appendix=="_apog1":
    apogeeflag = 1

    # what are the (hand-defined) components?
    comptag = dict()
    comptag[0] = ['disc','bar','knot']
    comptag[1] = ['disc','bar','knot']
    comptag[2] = ['disc','disc','bar']
    comptag[3] = ['-','disc','bar']

    # by component number, which component is [disc,bar,knot]?
    # assign the number of the cluster
    compnum = dict()
    compnum[0] = [0,1,2]
    compnum[1] = [0,1,2]
    compnum[2] = [0,1,-1]
    compnum[3] = [1,2,-1]

    # for ellipse tracing
    # by [bar,disc,knot]=[0,1,2], which component is which?
    # assign the number of the component
    complist = dict()
    complist[0] = [1,0,2]
    complist[1] = [1,0,2]
    complist[2] = [1,1,0]
    complist[3] = [-1,1,0]



# this should be repackaged
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


# where were these made?
# these are data points to make the angular momentum estimates: these can be plotted against lz for validation
modelradii = np.array([0.02,0.24,0.47,0.69,0.91,1.13,1.35,1.58,1.8,2.02,2.24,2.46,2.69,2.91,3.13,3.35,3.57,3.8,4.02,4.24,4.46,4.68,4.91,5.13,5.35,5.57,5.79,6.02,6.24,6.46,6.68,6.9,7.13,7.35,7.57,7.79,8.01,8.24,8.46,8.68,8.9,9.12,9.35,9.57,9.79,10.01,10.23,10.46,10.68,10.9,11.12,11.34,11.57,11.79,12.01,12.23,12.45,12.68,12.9,13.12,13.34,13.56,13.79,14.01,14.23,14.45,14.67,14.9,15.12,15.34,15.56,15.78,16.01,16.23,16.45,16.67,16.89,17.12,17.34,17.56,17.78,18.0,18.23,18.45,18.67,18.89,19.11,19.34,19.56,19.78,20.0,20.22,20.45,20.67,20.89,21.11,21.33,21.56,21.78,22.0])
barlz = np.array([0.01,1.52,5.55,12.1,21.17,32.76,46.87,63.5,82.65,104.32,128.51,155.21,184.44,216.19,250.46,287.25,326.56,368.39,412.73,459.6,508.99,560.9,615.33,672.27,731.74,793.73,858.24,925.26,994.81,1066.88,1141.47,1218.57,1298.2,1380.35,1465.01,1552.2,1641.91,1734.13,1828.88,1926.15,2025.93,2128.24,2233.06,2340.41,2450.28,2562.66,2677.57,2794.99,2914.94,3037.41,3162.39,3289.9,3419.92,3552.47,3687.53,3825.12,3965.22,4107.85,4252.99,4400.66,4550.84,4703.54,4858.77,5016.51,5176.78,5339.56,5504.87,5672.69,5843.03,6015.9,6191.28,6369.19,6549.61,6732.55,6918.02,7106.0,7296.5,7489.53,7685.07,7883.13,8083.71,8286.82,8492.44,8700.58,8911.24,9124.43,9340.13,9558.35,9779.09,10002.36,10228.14,10456.44,10687.26,10920.6,11156.46,11394.85,11635.75,11879.17,12125.11,12373.57])
massenclosedlz = np.array([0.02,22.02,78.4,138.56,194.93,248.5,300.17,349.9,397.97,444.17,489.49,533.72,578.7,624.52,672.7,722.94,776.78,832.63,891.37,951.85,1014.22,1078.71,1144.1,1211.36,1281.2,1352.71,1424.82,1498.96,1573.1,1648.77,1724.74,1801.5,1876.58,1952.4,2027.01,2102.42,2177.1,2252.43,2326.74,2400.58,2474.38,2548.83,2622.45,2694.79,2767.42,2840.21,2912.6,2983.37,3054.13,3126.27,3196.32,3266.51,3337.34,3407.65,3476.96,3546.65,3615.2,3685.63,3754.87,3823.76,3892.41,3960.97,4029.54,4099.28,4167.51,4236.03,4303.94,4370.83,4439.12,4507.8,4577.11,4645.08,4713.63,4780.89,4848.58,4915.95,4983.21,5050.05,5117.45,5186.44,5253.8,5321.42,5387.25,5452.63,5518.52,5585.23,5650.83,5716.79,5782.11,5848.28,5913.91,5977.2,6040.76,6104.13,6169.07,6234.38,6298.89,6362.73,6427.71,6491.23])
potentialmeasuredlz = np.array([0.24,15.72,54.65,106.62,162.29,217.05,272.69,330.3,385.12,433.17,476.67,519.94,561.94,598.29,627.95,656.25,690.39,734.31,785.7,839.32,890.91,939.58,986.11,1032.3,1079.91,1129.33,1181.09,1235.24,1289.75,1343.46,1396.86,1450.02,1502.46,1553.49,1604.04,1654.36,1704.45,1754.29,1803.89,1853.37,1902.8,1952.27,2001.74,2051.18,2100.68,2150.13,2199.53,2248.92,2298.23,2347.35,2396.26,2445.17,2493.98,2542.66,2591.35,2639.67,2688.28,2736.5,2784.52,2832.42,2879.86,2927.27,2974.35,3021.14,3067.84,3114.29,3160.31,3206.3,3252.34,3298.0,3343.3,3389.12,3434.67,3479.92,3525.01,3570.62,3616.0,3661.15,3706.09,3751.72,3797.22,3842.55,3887.72,3933.02,3978.71,4024.27,4069.68,4114.95,4160.31,4206.16,4251.89,4297.48,4342.93,4388.23,4433.86,4479.65,4525.3,4570.8,4616.16,4661.37])
