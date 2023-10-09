

inputdir   = 'data/apogee/'
appendix   = "" # only needed for model compatibility
modelappendix   = "_StarHorse" # only needed for model compatibility

modeltag = "APOGEE_all_feh8_SNRc_StarHorse"
datafile = "data/apogee/APOGEE_all_feh8_SNRc_StarHorse_tagged.txt"


modelname = "APOGEE"
mockanalysis = False

# what are the bins available to classify?
radii = [[0,1],[1,2],[2,3],[3,4],[4,5]]
binprefacs = [1.,1.,1.,1.,1.]


# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['knot','bar','disc']
comptag[1] = ['-','bar','disc']
comptag[2] = ['-','-','disc']
comptag[3] = ['disc','-','-']
comptag[4] = ['-','-','disc']

# for classifying
# by component number, which component is [disc,bar,knot]?
compnum = dict()
compnum[0] = [2,1,0]
compnum[1] = [2,1,-1]
compnum[2] = [2,-1,-1]
compnum[3] = [0,-1,-1]
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
