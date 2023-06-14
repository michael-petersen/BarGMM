
# what is the name of the mock?
inputdir   = 'data/apogee/'

# to recover the older classifications,
# https://github.com/michael-petersen/BarGMM/commit/1f50b03a9f6205e8b5d725458b4f83574bf109b6
#datafile = inputdir+"APOGEE_all_fehcut_reduceSNR_tagged.txt"
#modeltag = "APOGEE_all_fehcut_reduceSNR"

# https://github.com/michael-petersen/BarGMM/commit/dee2f4bcdd3e82e833a27616d74e3348786fe8fe
#modeltag = "APOGEE_all_feh6_SNRb"
#datafile = inputdir+"APOGEE_all_feh6_SNRb_tagged.txt"

modeltag = "APOGEE_all_feh6_SNRc"
datafile = inputdir+modeltag+"_tagged.txt"


modelname = "APOGEE"
mockanalysis = False

# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['disc','bar','knot']
comptag[1] = ['disc','bar','knot']
comptag[2] = ['disc','bar','knot']
comptag[3] = ['disc','bar','-']

# for classifying
# by component number, which component is [disc,bar,knot]?
compnum = dict()
compnum[0] = [0,1,2]
compnum[1] = [0,1,2]
compnum[2] = [0,1,2]
compnum[3] = [0,1,-1]

# for ellipse tracing
# by [bar,disc,knot], which component is which index?
# (the first component is which component by the tags above?)
complist = dict()
complist[0] = [1,0,2]
complist[1] = [1,0,2]
complist[2] = [1,0,2]
complist[3] = [1,0,-1]


# half-integer bins
comptag[5] = ['disc','bar','knot']
comptag[15] = ['disc','bar','knot']
comptag[25] = ['disc','bar','-']
comptag[35] = ['-','-','disc']

compnum[5] = [0,1,2]
compnum[15] = [0,1,2]
compnum[25] = [0,1,-1]
compnum[35] = [2,-1,-1]

complist[5] = [1,0,2]
complist[15] = [1,0,2]
complist[25] = [1,0,-1]
complist[35] = [-1,-1,0]
