
import numpy as np

inputdir     = 'data/barmock/'
datafile     = inputdir+"ModelBarYoungMock10000smallerror.txt"
modeltag     = "ModelBarYoungMock10000smallerror"
appendix     = "_apog1"  # if this is for the SDSS-coverage sample

modelname    = 'BarMock'
mockanalysis = True


apogeeflag = 1

radii      = [[0,1],[1,2],[15,25],[2,3],[25,35],[3,4]]
binprefacs = [ 1.,   1.,   0.1,    1.,   0.1,   1.]

# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['disc','bar','knot']
comptag[1] = ['disc','disc','bar']
comptag[2] = ['disc','-','-'] # in this bin, the bar has been split into two components for some reason.
comptag[3] = ['-','disc','bar']
comptag[15] = ['disc','disc','bar'] # in this bin, the bar has been split into two components for some reason.
comptag[25] = ['-','disc','bar'] # in this bin, the bar has been split into two components for some reason.
comptag[35] = ['-','-','disc']

# by component number, which component is [disc,bar,knot]?
# assign the number of the cluster
compnum = dict()
compnum[0] = [0,1,2]
compnum[1] = [0,2,-1]
compnum[2] = [0,-1,-1]
compnum[3] = [1,2,-1]
compnum[15] = [0,2,-1]
compnum[25] = [1,2,-1]
compnum[35] = [2,-1,-1]
