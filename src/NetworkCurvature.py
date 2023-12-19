from typing import Dict
import matplotlib.pyplot as plt
import sys
import numpy as np
import networkx as nx
import ot
import tqdm 
import time

def compute_node_curvatures(
	G:nx.Graph,
	weight:str = 'weight',
	edge_curv:str = 'ORC',
	node_curv:str = 'node_curvature',
	norm_node_curv:str  = 'normalized_curvature'
	)->None:
	
	for node in G.nodes():
		normalized_curvature = 0 
		raw_curvature = 0 
		for x in G.neighbors(node):
			term_weight = G[node][x][weight]
			sum_term = G[node][x][edge_curv]
			raw_curvature += sum_term
			normalized_curvature += sum_term*term_weight
		nx.set_node_attributes(G,{node:raw_curvature}, node_curv)
		nx.set_node_attributes(G,{node:normalized_curvature},norm_node_curv)

def assign_densities(
	G:nx.Graph,
	alpha:float = 0.5,
	weight:str = 'weight',
	measure_name:str = 'density'
	) -> None:
	
	assert 0<=alpha<=1, 'alpha must be between 0 and 1'
	attrs = nx.get_node_attributes(G,'Gene')
	for node in G.nodes():
		density = np.zeros(len(G.nodes()))
		density[node] = alpha
		node_degree = G.degree(node,weight = weight)
		# print("NODE:\t{n}\tDEGREE:\t{d}".format(n=node, d=node_degree))

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

	

	path_lengths = dict(nx.all_pairs_dijkstra_path_length(G,weight=weight))
	for node1 in path_lengths.keys():
		node1Lengths = path_lengths[node1]
		for node2 in node1Lengths.keys():
			D[node1,node2] = np.round(node1Lengths[node2],5)
			# if D[node1,node2]==0 and node1!=node2:
			# 	print(node1)
			# 	print(node2)
			# 	print(node1Lengths[node2])
			# 	print(node1Lengths)
			# 	sys.exit(1)
							 #rounding to make sure D is symmetric
			
	# # slow b/c redundant
	# for n1 in G.nodes():
	# 	for n2 in G.nodes():
	# 		D[n1,n2] = np.round(nx.dijkstra_path_length(G,n1,n2,weight = weight),7)

	if not (D==D.T).all():
		print('symmetry error')
		print(D==D.T)
		issues = np.where(D!=D.T)
		print(D[issues[0][0],issues[1][0]])
		print(D[issues[1][0],issues[0][0]])
		sys.exit(1)

	return np.ascontiguousarray(D)

	




def compute_OR_curvature(
	G:nx.Graph,
	D:np.ndarray,
	density:str = 'density',
	curvature_name:str =  'ORC',
	sinkhorn: bool = False,
	epsilon: float = 0.01
	) -> None:
	

	assert epsilon>0, 'epsilon must be positive'

	
	
	


	for edge in G.edges():
		u,v = edge[0],edge[1]
		m_u = G.nodes[u][density]
		m_v = G.nodes[v][density]


		if sinkhorn:
			
			W = ot.sinkhorn2(a= m_u, b= m_v,M= D,reg = epsilon,warn = False)
			
			
		else:
			W = ot.emd2(a= m_u, b= m_v, M =D, numItermax=1000000)

		kappa =  1- (W/D[u,v])
		G[u][v][curvature_name] = np.round(kappa,2)

		

		

	
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


	


	G = nx.barbell_graph(4,1)
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


	
