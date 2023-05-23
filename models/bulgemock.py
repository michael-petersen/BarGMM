

inputdir  = 'data/bulgemock/'
datafile  = inputdir+"ModelBarYoungBulgeMock10000.txt"
modeltag  = "ModelBarYoungBulgeMock10000"
modelname = 'BulgeMock'

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
