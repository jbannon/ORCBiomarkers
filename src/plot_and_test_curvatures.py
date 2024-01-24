import sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np 
import NetworkCurvature as nc
import os


# n_iter = 200
# count = 0 
rng = np.random.seed(1234)
# for i in range(n_iter):
# 	G = nx.fast_gnp_random_graph(10,0.3,rng)
# 	orc = nc.OllivierRicciCurvature(G)
# 	orc.compute_edge_curvatures()
# 	G_ = orc.G.copy()
# 	my_dict = {}
# 	for e in G_.edges(data=True):
# 		key = (e[0],e[1])
# 		value = e[2]['ricci_curvature']
# 		my_dict[key] = value

# 	orc = OllivierRicci(G)
# 	orc.compute_ricci_curvature()
# 	G_ = orc.G.copy()
# 	their_dict = {}

# 	for e in G_.edges(data=True):
# 		key = (e[0],e[1])
# 		value = e[2]['ricciCurvature']
# 		their_dict[key] = value

# 	if their_dict == my_dict:
# 		count+=1


# print(count)


def plot_curvature_hist(graph_label, curvatures):
	plt.hist(curvatures,bins=20)
	plt.xlabel('Ricci curvature')
	plt.ylabel('Count')
	plt.title("{g} Graph Edge Curvature Histogram".format(g=graph_label))
	os.makedirs("../figs/examples/",exist_ok = True)
	plt.savefig("../figs/examples/{g}_hist.png".format(g=graph_label))
	plt.close()


def plot_graph(
	G,
	graph_label,
	curvatures
	):
	colors=curvatures
	cmap=plt.cm.viridis
	vmin = min(curvatures)
	vmax = max(curvatures)
	if graph_label == 'Grid':
		pos = {(x,y):(y,-x) for x,y in G.nodes()}
	# elif graph_label == 'Tree':
	# 	pos = nx.nx_agraph.graphviz_layout(G)
	else:
		pos=nx.circular_layout(G)
		pos = nx.spring_layout(G,k=5,pos=pos,seed=12345)
	fig, ax = plt.subplots()
	nx.draw(G,pos, node_color='black',edge_color=curvatures, width=4, edge_cmap=cmap,with_labels=False, vmin=vmin, vmax=vmax)
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = vmin, vmax=vmax))
	cbar = plt.colorbar(sm)#, fraction=0.5, pad=0)
	cbar.ax.set_ylabel('Edge Curvature',labelpad = 20,rotation = 270,fontsize=16)
	plt.savefig("../figs/examples/{g}_graph.png".format(g=graph_label))
	plt.close()
	
	
for graph_label, G in zip(['Barbell','Grid','Random','Complete','Tree'],
	[nx.barbell_graph(5,0),nx.grid_2d_graph(8, 8),
	nx.fast_gnp_random_graph(15,0.5,rng),nx.complete_graph(9),
	nx.random_tree(6)]):
	orc = nc.OllivierRicciCurvature(G)
	orc.compute_edge_curvatures()
	G = orc.G.copy()

	edges, curvatures = zip(*nx.get_edge_attributes(G,"ricci_curvature").items())
	plot_curvature_hist(graph_label, curvatures)
	plot_graph(G,graph_label,curvatures)







