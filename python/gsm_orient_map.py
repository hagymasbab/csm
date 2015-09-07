from igraph import *
from numpy import *
from scipy.io import loadmat

filenames = ["cutoff01","cutoff02","cutoff03","cutoff04"]

for fn in filenames:
	print "Inputfile " + fn
	
	corrmat = loadmat("c_adj_%s.mat" % fn)["CC"]
	A = ndarray.tolist(corrmat) 
	g = Graph.Weighted_Adjacency(A,mode=ADJ_UNDIRECTED,loops=False)
	#g = g.simplify()

	ors = loadmat('../matlab/bin/orientpref_gsm-learned_248.mat')['orientpref'][0]
	colors = []
	centers = [0,22.5,45,67.5,90,112.5,135,157.5,180]
	colornames = ['cyan','green','lightgreen','yellow','red','pink','purple','blue','cyan']
	for i in range(0,g.vcount()):
		for j in range(len(centers)-1):
			if ors[i] >= centers[j] and ors[i] < centers[j+1]:
				if ors[i] - centers[j] >= centers[j+1] - ors[i]: colors.append(colornames[j])
				else: colors.append(colornames[j+1])
		#if ors[i] > 28 and ors[i] < 65: colors.append("red")
		#elif ors[i] >= 65 and ors[i] < 85: colors.append("blue")
		#elif ors[i] >= 85 and ors[i] < 110: colors.append("green")
		#else: colors.append("black")
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
	osi_thresholds = [0.1,0.2]	
	for ot in osi_thresholds:
		print "..OSI threshold %.1f" % ot
		to_delete = []
		for i in range(0,g.vcount()):
			if osi[i] < ot: to_delete.append(i)
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
		layouts["mdsw"] = g.layout_mds(dist=distmat)
		layouts["mds"] = g.layout_mds()
		layouts["graphopt"] = g.layout_graphopt()
		
		for lk in layouts.keys():
			print "....Layout " + lk
			outname = "orient_map_%s_osith%.1f_%s.svg" % (fn,ot,lk);
			g.write_svg(outname,layout=layouts[lk],labels=None,vertex_size=10,colors="color",edge_colors=["grey"]*g.ecount())

