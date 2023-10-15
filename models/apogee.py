
inputdir   = 'data/apogee/'
appendix   = "" # only needed for model compatibility
modelappendix = ""

# version using astroNN distances
modeltag = "APOGEE_all_feh6_SNRc"
datafile = inputdir+modeltag+"_tagged.txt"

# version using Bailer-Jones distances
modeltag = "APOGEE_all_feh8_SNRc"
datafile = inputdir+modeltag+"_tagged.txt"


modelname = "APOGEE"
mockanalysis = False

# what are the bins available to classify?
radii = [[0,1],[5,15],[1,2],[15,25],[2,3],[25,35],[3,4],[35,45],[4,5]]
binprefacs = [1.,0.1,1.,0.1,1.,0.1,1.,0.1,1.0]


# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['disc','bar','knot']
comptag[1] = ['disc','bar','knot']
comptag[2] = ['disc','bar','knot']
comptag[3] = ['disc','-','bar']
comptag[4] = ['-','-','disc']

# for classifying
# by component number, which component is [disc,bar,knot]?
compnum = dict()
compnum[0] = [0,1,2]
compnum[1] = [0,1,2]
compnum[2] = [0,1,2]
compnum[3] = [0,2,-1]
compnum[4] = [2,-1,-1]

# for ellipse tracing
# by [bar,disc,knot], which component is which index?
# (the first component is which component by the tags above?)
complist = dict()
complist[0] = [1,0,2]
complist[1] = [1,0,2]
complist[2] = [1,0,2]
complist[3] = [1,-1,0]
complist[4] = [-1,-1,0]


# half-integer bins
comptag[5] = ['disc','bar','knot']
comptag[15] = ['disc','bar','knot']
comptag[25] = ['-','disc','bar']
comptag[35] = ['-','-','disc']

compnum[5] = [0,1,2]
compnum[15] = [0,1,2]
compnum[25] = [1,2,-1]
compnum[35] = [2,-1,-1]

complist[5] = [1,0,2]
complist[15] = [1,0,2]
complist[25] = [-1,1,-0]
complist[35] = [-1,-1,0]
