import matplotlib.pyplot as plt
import argparse
from typing import Dict, Union, List
import sys
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import numpy as np 
import networkx as nx
import ot
import yaml 
import utils
# import statsmodels.stats.multitest.multipletests
import pickle

import statsmodels.stats.multitest as mt



def main(
	Config:Dict
	) -> None:
	
	drug:str = 'Atezo'
	tissue:str = 'BLCA'
	# alphas:List[float] 
	alpha:float = 0.5
	# num_iters:int
	# expr_dir:str
	# geneset_dir:str
	# do_one_hop:bool
	# norm: str
	rng_seed: int = 12345
	# mode: str
	# fdr_threshold: float
	# count_cutoff: int
	# edge_type: str ME, C  (mass action/ correlation)
	# curvature_type: str OR
	# connectivity: str = 'sparse' 
	# weighting: str = 'unweighted'
	network_type:str = 'ME'
	np.random.seed(rng_seed)
	# drug, tissue, alphas, num_iters, expr_dir, geneset_dir, do_one_hop, norm =\
	# 	utils.UNPACK_PARAMETERS(config['EXPERIMENT_PARMS']) 
	
	# fdr_threshold, count_cutoff = utils.UNPACK_PARAMETERS(config['DE_PARAMS'])



	# utils.make_file_name(['a','b','c'])
	# expr_file = "{data_directory}/{drug}/{tissue}/expression.csv".format()
	expr = pd.read_csv("../data/expression/{d}/{t}/expression.csv".format(d=drug,t=tissue))
	resp = pd.read_csv("../data/expression/{d}/{t}/response.csv".format(d=drug,t=tissue))
	DE = pd.read_csv("../data/genesets/{d}_{t}_DE.csv".format(d=drug,t=tissue),index_col = 0)
	
	netfile = "../data/networks/{connectivity}/{weighting}.pickle".format(connectivity = 'sparse',weighting = 'weighted')
	
	with open(netfile,"rb") as istream:
		PPI_Graph = pickle.load(istream)
	
	do_one_hop = False

	DE = DE[DE['Thresh.Value']==0.05]
	DE = DE[DE['Count']>=150]
	
	# DE_genes = list(DE['Gene'].values)
	idleness = 0.5
	
	gene_list = list(DE['Gene'].values)[:100]

	common_genes = list(set(gene_list).intersection(set(expr.columns[1:])))
	
	# issues = list(set(expression_data.columns[1:]).difference(set(gene_list)))
	print("Gene list length: {L}".format(L =len(gene_list)))
	print("Common Genes length: {L}".format(L =len(common_genes)))


	if do_one_hop:
		seeds = [x for x in common_genes if x in PPI_Graph.nodes()]
		one_hop = [x for x in seeds]
		for n in seeds:
			nbrs = PPI_Graph.neighbors(n)
			one_hop.extend([x for x in nbrs])
		common_genes = one_hop
	

	expr = expr[['Run_ID']+common_genes]
	
	responder_IDs = resp[resp['Response']==1]['Run_ID'].tolist()
	num_samples = 10
	LCC_Graph = utils.harmonize_graph_and_geneset(PPI_Graph,common_genes)
	print(len(LCC_Graph))
	responder_IDs = resp[resp['Response']==1]['Run_ID'].tolist()
	responder_samples = np.random.choice(responder_IDs, num_samples, replace = False)
	responder_measurements = expr[expr['Run_ID'].isin(responder_samples)]
	print(responder_measurements)
	pvals = []
	# revisit graph LCC creation to preserve order. 
	# dumb solution is to make a new graph
	# and add elements in to it with node ids and weights. 
	# need a way to assign weights to neighbors. 
	# then compute OR curvature. 
	# focus on wavelet moments and magnetic laplacian wavelet
	# GL = list(LCC_Graph.nodes)
	
	# 	v1 = responder_measurements[g1].values
	# 	v2 = responder_measurements[g2].values
	# 	weight, pval = stats.spearmanr(v1,v2)
	# 	pvals.append(pval)

	# 	LCC_Graph[g1][g2]['weight'] = 0.5*(1+ weight)
	# 	LCC_Graph[g1][g2]['ORC'] = None
	
	G = nx.Graph()
	GL = list(LCC_Graph.nodes)
	node_to_idx  = {}
	for i in range(len(GL)):
		G.add_node(i,gene_name = GL[i],total_curvature = None, density = np.zeros(len(GL)))
		node_to_idx[GL[i]] = i

	for edge in LCC_Graph.edges():
		g1, g2 = edge[0],edge[1]
		v1 = responder_measurements[g1].values
		v2 = responder_measurements[g2].values
		weight, pval = stats.spearmanr(v1,v2)
		weight = 0.5*(1+weight)
		pvals.append(pval)
		G.add_edge(node_to_idx[g1],node_to_idx[g2],weight=weight, ORC = None)
	
	# print(G)
	# nx.draw(G)
	# plt.show()

	# for n in G.nodes(data=True):
	# 	dx = G.degree(n,weight='weight')
	# 	print(dx)
	# 	# print(n[1]['density'])
	# 	print(n)
	


	# node_to_idx = {}
	# idx_to_node = {}
	
	# for i in range(len(GL)):
	# 	node_to_idx[GL[i]] = i
	# 	idx_to_node[i] = GL[i]

	densities = {i:np.zeros(len(GL)) for i in range(len(GL))}


	
	

	
	# if do_one_hop:
	# 	pass
	# 	# collect the neighbors then reduce to the LCC

	
	
	# print(nx.is_connected(G_0))
	# assign_densities(G_0)



if __name__ == '__main__':
	# parser = argparse.ArgumentParser()
	# parser.add_argument("-config",help="The config file for these experiments")
	# args = parser.parse_args()
	
	# with open(args.config) as file:
	# 	config = yaml.safe_load(file)

	# main(config)
	main({})
	