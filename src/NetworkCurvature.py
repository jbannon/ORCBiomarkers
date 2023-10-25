from typing import Dict
import matplotlib.pyplot as plt
import sys
import numpy as np
import networkx as nx
import ot




def assign_densities(
	G:nx.Graph,
	alpha:float = 0.5,
	weight:str = 'weight',
	measure_name:str = 'density'
	) -> None:
	
	assert 0<=alpha<=1, 'alpha must be between 0 and 1'
	
	for node in G.nodes():
		density = np.zeros(len(G.nodes()))
		density[node] = alpha
		node_degree = G.degree(node,weight = weight)
		
		for x in G.neighbors(node):
			density[x] = (1-alpha)*( (G[node][x][weight])/node_degree)

		nx.set_node_attributes(G,{node:np.ascontiguousarray(density)},measure_name)


	

def make_APSP_Matrix(
	G:nx.Graph,
	weighted:bool = True,
	weight:str = 'weight'
	) -> np.ndarray:
	

	N = len(G.nodes)
	D = np.zeros((N,N))
	
	if not weighted:
		weight = None
	# slow b/c redundant
	for n1 in G.nodes():
		for n2 in G.nodes():
			D[n1,n2] = np.round(nx.dijkstra_path_length(G,n1,n2,weight = weight),4)

	if not (D==D.T).all():
		print('symmetry error')

	return np.ascontiguousarray(D)

	




def compute_OR_curvature(
	G:nx.Graph,
	D:np.ndarray,
	density:str = 'vertex_measure',
	weight:str = 'weight',
	node_to_idx:Dict = None,
	curvature_field:str =  'ORC',
	sinkhorn: bool = False,
	epsilon: float = 0.01
	) -> None:
	
	assert epsilon>0, 'epsilon must be positive'

	for edge in G.edges():
		u,v = edge[0],edge[1]
		m_u = G.nodes[u][density]
		m_v = G.nodes[v][density]
		if sinkhorn:
			W = ot.sinkhorn(a= m_u, b= m_v,M= D,reg = epsilon)
		else:
			
			W = ot.emd2(a= m_u, b= m_v, M =D)
		G[u][v][curvature_field] = np.round(1- (W/D[u,v]),3)
		

		

	
if __name__ == '__main__':
	rng = np.random.default_rng(seed=12345)
	G = nx.fast_gnp_random_graph(10,0.3,rng)
	np.random.seed(12345)
	for edge in G.edges():
		G[edge[0]][edge[1]]['weight'] = np.random.uniform(1,2)
		G[edge[0]][edge[1]]['ORC'] = None
		

	D = make_APSP_Matrix(G)
	assign_densities(G)
	compute_OR_curvature(G,D)

	pos = nx.spring_layout(G)
	nx.draw(G,pos)
	labels = nx.get_edge_attributes(G,'ORC')
	nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
	plt.show()


	


	G = nx.barbell_graph(4,0)
	for edge in G.edges():
		G[edge[0]][edge[1]]['weight'] = np.random.uniform(1,2)
		G[edge[0]][edge[1]]['ORC'] = None
	D = make_APSP_Matrix(G)
	assign_densities(G)
	compute_OR_curvature(G,D)
	edges,weights = zip(*nx.get_edge_attributes(G,'ORC').items())


	pos=nx.circular_layout(G)
	pos = nx.spring_layout(G,k=5,pos=pos,seed=12345)
	nx.draw(G, pos,node_color='b', edgelist=edges, edge_color=weights,width = 2.5,with_labels=True)
	labels = nx.get_edge_attributes(G,'ORC')

	nx.draw_networkx_edge_labels(G,pos,edge_labels=labels,bbox=dict(boxstyle = 'round',facecolor='white',edgecolor='black'),rotate=True)
	plt.legend()
	plt.show()


	
