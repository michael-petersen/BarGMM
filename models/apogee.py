
# what is the name of the mock?
inputdir   = 'data/apogee/'
datafile = inputdir+"APOGEE_all_fehcut_reduceSNR_tagged.txt"
modeltag = "APOGEE_all_fehcut_reduceSNR"
modelname = "APOGEE"
mockanalysis = False

# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['bar','disc','knot']
comptag[1] = ['disc','bar','knot']
comptag[2] = ['disc','-','bar']
comptag[3] = ['-','-','disc']
comptag[4] = ['-','-','disc']

# by component number, which component is [disc,bar,knot]?
compnum = dict()
compnum[0] = [1,0,2]
compnum[1] = [0,1,2]
compnum[2] = [0,2,-1]
compnum[3] = [2,-1,-1]
compnum[4] = [2,-1,-1]

# for ellipse tracing
# by [bar,disc,knot], which component is which index?
complist = dict()
complist[0] = [0,1,2]
complist[1] = [1,0,2]
complist[2] = [1,-1,0]
complist[3] = [-1,-1,1]
complist[4] = [-1,-1,1]


# half-integer bins
comptag[5] = ['disc','bar','knot']
comptag[15] = ['-','-','bar']
comptag[25] = ['-','disc','bar']

compnum[5] = [0,1,2]
compnum[15] = [-1,2,-1]
compnum[25] = [1,2,-1]

complist[5] = [1,0,2]
complist[15] = [2,-1,-1]
complist[25] = [2,1,-1]
