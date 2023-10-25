from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split, GridSearchCV
from collections import defaultdict
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
import os
import yaml 
import utils
# import statsmodels.stats.multitest.multipletests
import pickle
import NetworkCurvature as nc
import statsmodels.stats.multitest as mt
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.linear_model import LogisticRegression
import tqdm








def main(
	Config:Dict
	) -> None:
	
	drug:str 
	tissue:str
	dataDir:str
	genesetDir:str
	doOneHop:bool
	rngSeed: int 
	fdr_threshold: float
	count_cutoff: int
	graphTops:List[str]
	weighting: str
	
	
	np.random.seed(rng_seed)
	
	

	drug, tissue, alphas, rngSeed,dataDir, genesetDir, doOneHop = \
		utils.UNPACK_PARAMETERS(config['EXPERIMENT_PARMS']) 
	
	# drug = 'Atezo'
	# tissue = 'BLCA'
	
	exprFile = utils.make_file_path(dataDir,[drug, tissue],'expression','.csv')
	respFile = utils.make_file_path(dataDir,[drug,tissue],'response','.csv')


	resp = pd.read_csv(respFile)
	expr = pd.read_csv(exprFile)
	
	
	# fdr_threshold, count_cutoff = utils.UNPACK_PARAMETERS(config['DE_PARAMS'])



	# utils.make_file_name(['a','b','c'])
	# expr_file = "{data_directory}/{drug}/{tissue}/expression.csv".format()
	
	
	
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

	common_genes = [x for x in gene_list if x in expr.columns[1:]]

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
	
	expr = expr.merge(resp, on = 'Run_ID')
	
	
	LCC_Graph = utils.harmonize_graph_and_geneset(PPI_Graph,common_genes)

	
	
	
	
	
	# print(GL)
	gene_to_idx  = {} # cahnge to 'gene to idx'
	idx_to_gene = {}
	for idx, gene in enumerate(LCC_Graph.nodes):
		gene_to_idx[gene] = idx
		idx_to_gene[idx] = gene
	
	LCC_Graph = nx.relabel_nodes(LCC_Graph,gene_to_idx)
	nx.set_node_attributes(LCC_Graph,idx_to_gene, 'Gene')
	

	edge_curvatures = defaultdict(list)
	node_curvatures = defaultdict(list)
	
	for idx, row in expr.iterrows():
		G = LCC_Graph.copy()
		r = row['Response']
		
		total_edge_weight = 0
		for edge in G.edges():
			v1, v2 = edge[0], edge[1]
			g1, g2 = idx_to_gene[v1], idx_to_gene[v2]
			w = np.round(np.log2(row[g2]+1.00001)*np.log2(row[g1]+1),2)
			total_edge_weight += w
			G[v1][v2]['weight'] = w

		for edge in G.edges():
			v1, v2 = edge[0],edge[1]
			G[v1][v2]['weight'] = 10*G[v1][v2]['weight']/total_edge_weight
			
		nc.assign_densities(G)
		D = nc.make_APSP_Matrix(G)
		nc.compute_OR_curvature(G,D,'density')
		for edge in G.edges(data=True):
			
			e = "{a}/{b}".format(a= idx_to_gene[edge[0]],b = idx_to_gene[edge[1]])
			edge_curvatures['edge'].append(e)
			edge_curvatures['curvature'].append(edge[2]['ORC'])
			edge_curvatures['Response'].append(r)
		

	df = pd.DataFrame(edge_curvatures)


	keep_edges = []
	pvals = []
	for e in pd.unique(df['edge']):
		tempR = df[df['edge']==e]
		tempR = tempR[tempR['Response']==1]
		tempNR = df[df['edge']==e]
		tempNR = tempNR[tempNR['Response']==0]
		# print(tempR)
		# print(tempNR)
		u, p = mannwhitneyu(tempR['curvature'].values,tempNR['curvature'].values,alternative = 'two-sided')
		keep_edges.append(e)
		pvals.append(p)

	res = mt.multipletests(pvals,method = 'fdr_bh')
	adjust_pvals = res[1]
	# print()
	keep_edges = [keep_edges[i] for i in np.where(adjust_pvals<=0.1)[0]]
	print(keep_edges)

	df = df[df['edge'].isin(keep_edges)]

	# sns.boxplot(x = 'edge', y = 'curvature',hue = 'Response',data = df)
	# plt.show()
	# print(keep_edges)
	#['KIF2C/CDC25C', 'KIF2C/CENPF', 'MAD2L1/SPC25', 'MCM6/RFC3']
	
	
	# v1 = responder_measurements[g1].values
	# v2 = responder_measurements[g2].values
		
	
	


	
	

	
	



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)
	# main({})
	