

inputdir  = 'data/barmock/'
datafile  = inputdir+"ModelBarYoungMock10000.txt"
modeltag  = "ModelBarYoungMock10000"
modelname = 'BarMock'
mockanalysis = True

# what are the (hand-defined) components?
comptag = dict()
comptag[0] = ['disc','bar','bar']
comptag[1] = ['disc','disc','bar']
comptag[2] = ['disc','-','bar']
comptag[3] = ['-','-','disc']
comptag[4] = ['-','-','disc']

# by component number, which component is [disc,bar,knot]?
# assign the number of the cluster
compnum = dict()
compnum[0] = [0,2,-1]
compnum[1] = [1,2,-1]
compnum[2] = [0,2,-1]
compnum[3] = [2,-1,-1]
compnum[4] = [2,-1,-1]

# for ellipse tracing
# by [bar,disc,knot]=[0,1,2], which component is which?
# assign the number of the component
complist = dict()
complist[0] = [1,0,0]
complist[1] = [1,1,0]
complist[2] = [1,-1,0]
complist[3] = [-1,-1,1]
complist[4] = [-1,-1,1]
