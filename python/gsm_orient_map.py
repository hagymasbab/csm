from igraph import *
from numpy import *
from scipy.io import loadmat

filenames = ["cutoff02","cutoff03","cutoff04"]

for fn in filenames:
	print fn
	
	corrmat = loadmat("c_adj_%s.mat" % fn)["CC"]
	A = ndarray.tolist(corrmat) 
	g = Graph.Weighted_Adjacency(A,mode=ADJ_UNDIRECTED)
	g = g.simplify()

	ors = loadmat('../matlab/bin/orientpref_gsm-learned_248.mat')['orientpref'][0]
	colors = []
	for i in range(0,g.vcount()):
		if ors[i] > 28 and ors[i] < 65: colors.append("red")
		elif ors[i] >= 65 and ors[i] < 85: colors.append("blue")
		elif ors[i] >= 85 and ors[i] < 110: colors.append("green")
		else: colors.append("black")
		#col1d = (ors[i]/180)*255
		#colors.append("#%02X%02X%02X" % (col1d,0,255-col1d))
	g.vs["color"] = colors
	
	distances = []
	for e in g.es:
		corrval = corrmat[e.source,e.target]
		if corrval == 0: distances.append(1)
		else: distances.append(1/corrval)
	g.es["distance"] = distances

	osi = loadmat('../matlab/bin/osi_gsm-learned_248.mat')['OSI'][0]
	to_delete = []
	for i in range(0,g.vcount()):
		if osi[i] < 0.2: to_delete.append(i)
	g.delete_vertices(to_delete)
	
	g = g.clusters().giant()
	
	distmat = ones((g.vcount(),g.vcount()))
	for e in g.es:
		distmat[e.source,e.target] = e["distance"]
		distmat[e.target,e.source] = e["distance"]
	
	layouts = {}
	layouts["fr"] = g.layout_fruchterman_reingold(weights="weight",area=16*g.vcount()**2)
	#layouts["drl"] = g.layout_drl(weights="weight")
	#lay = g.layout_kamada_kawai(weights="weight")
	#lay = g.layout_sugiyama(weights="weight")
	#layouts["mdsw"] = g.layout_mds(dist=distmat)
	#layouts["mds"] = g.layout_mds()
	layouts["graphopt"] = g.layout_graphopt()
	
	for lk in layouts.keys():
		print lk
		outname = "orient_map_%s_%s.svg" % (fn,lk);
		g.write_svg(outname,layout=layouts[lk],labels=None,vertex_size=10,colors="color",edge_colors=["grey"]*g.ecount())

