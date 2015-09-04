from igraph import *
from numpy import *
from scipy.io import loadmat

A = ndarray.tolist(loadmat("c_adj.mat")["CC"]) 
#A = ndarray.tolist(loadmat("../matlab/bin/c_gsm-learned_248.mat")["C"]) 
g = Graph.Weighted_Adjacency(A,mode=ADJ_UNDIRECTED,loops=False)
ors = loadmat('../matlab/bin/orientpref_gsm-learned_248.mat')['orientpref'][0]
colors = []
for i in range(0,g.vcount()):
	col1d = (ors[i]/180)*255
	colors.append("#%02X%02X%02X" % (col1d,0,255-col1d))
g.vs["color"] = colors
g = g.clusters().giant()
lay = g.layout_fruchterman_reingold(weights="weight")
g.write_svg('orient_map.svg',layout=lay,labels=None,vertex_size=5,colors="color")

