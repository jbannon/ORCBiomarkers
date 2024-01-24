import tqdm
import sys
import pandas as pd
import numpy as np 
import networkx as nx
from typing import List, Dict
import statsmodels.stats.multitest as mt
from scipy.stats import mannwhitneyu, wilcoxon
from collections import defaultdict
import NetworkCurvature as nc






class CurvatureFeatureSelector:
	def __init__(
		self,
		base_graph:nx.Graph,
		X:np.ndarray,
		y:np.ndarray,
		gene_to_col:Dict[int,str],
		gene_to_node_num:Dict[str,int],
		node_num_to_gene:Dict[int,str],
		alpha:float = 0.5,
		pvalue_thresh:float = 0.1,
		weight_field:str = 'weight',
		path_method:str = "pairwise",
		curvature_field:str = "ricci_curvature",
		node_field:str = "node_curvature",
		norm_node_field:str = "node_curvature_normalized",
		measure_name:str = 'density',
		min_distance:float = 1E-5,
		min_degree:float = 1E-5,
		sinkhorn_thresh:int = 2000,
		epsilon:float = 1.0
		) -> None:
		
		assert X.shape[0] == y.shape[0], "X and Y must have equal first dim"

		self.G = base_graph.copy()
		
		self.X = X
		self.y = y 
		self.gene_to_col = gene_to_col
		self.gene_to_node_num = gene_to_node_num
		self.node_num_to_gene = node_num_to_gene
		self.alpha = alpha
		self.pvalue_thresh = pvalue_thresh
		self.weight_field = weight_field
		self.path_method = path_method
		self.curvature_field = curvature_field
		self.node_field = node_field
		self.norm_node_field = norm_node_field,
		self.measure_name = measure_name
		self.min_distance = min_distance
		self.min_degree = min_degree
		self.sinkhorn_thresh = sinkhorn_thresh
		self.epsilon = epsilon

		self.edge_curvatures = defaultdict(list)
		self.node_curvatures = defaultdict(list)			
		self.raw_node_curvatures = defaultdict(list)

		
		



	def compute_curvatures(
		self
		) -> None:

		
		for i in tqdm.tqdm(range(self.X.shape[0])):
			sample_response = self.y[i]
			
			for edge in self.G.edges():
				node1, node2 = edge[0], edge[1]
				gene1, gene2 = self.node_num_to_gene[node1], self.node_num_to_gene[node1]
				col1, col2 = self.gene_to_col[gene1],self.gene_to_col[gene2]
				weight = np.round(np.log2(self.X[i,col1]+1)*np.log2(self.X[i,col2]+1),5)
				self.G[node1][node2][self.weight_field] = weight

			orc = nc.OllivierRicciCurvature(
				G= self.G,
				alpha = self.alpha, 
				weight_field = self.weight_field,
				path_method = self.path_method,
				curvature_field = self.curvature_field,
				node_field = self.node_field,
				norm_node_field = self.norm_node_field,
				measure_name = self.measure_name,
				min_distance = self.min_distance,
				min_degree = self.min_degree,
				sinkhorn_thresh = self.sinkhorn_thresh,
				epsilon = self.epsilon
				)
			
			orc.compute_curvatures()
			G = orc.G.copy()
		
		

			for edge in G.edges(data=True):
				self.edge_curvatures['Gene1'].append(self.node_num_to_gene[edge[0]])
				self.edge_curvatures['Gene2'].append(self.node_num_to_gene[edge[1]])
				edge_string= "{a};{b}".format(a=self.node_num_to_gene[edge[0]] ,b = self.node_num_to_gene[edge[1]])				
				self.edge_curvatures['Edge'].append(edge_string)
				self.edge_curvatures['Curvature'].append(edge[2][self.curvature_field])
				self.edge_curvatures['Response'].append(sample_response)
			
			for node in G.nodes(data= True):
				self.node_curvatures['Gene'].append(node[1]['gene'])
				self.node_curvatures['Raw'].append(node[1][self.node_field])
				self.node_curvatures['Normalized'].append(node[1][self.norm_node_field])
				self.node_curvatures['Response'].append(sample_response)

		
		

	def get_features(
		self,
		feature_type:str
		):

		if feature_type == "edge":
			genes, idxs = self.compute_edge_pvalues()
		elif feature_type == "node":
			genes, idxs = self.compute_node_pvalues()
		elif feature_type == "node-norm":
			genes, idxs = self.compute_normalized_node_pvalues()

		return genes, idxs
	

	def compute_edge_pvalues(self):
		
		edge_data = pd.DataFrame(self.edge_curvatures)
		edge_res = defaultdict(list)
		edge_pvals = []

		for e in pd.unique(edge_data['Edge']):
			tempR = edge_data[edge_data['Edge'] == e]
			tempR = tempR[tempR['Response'] == 1]
			tempNR = edge_data[edge_data['Edge'] == e]
			tempNR = tempNR[tempNR['Response'] == 0 ]
			u, p = mannwhitneyu(tempR['Curvature'].values,tempNR['Curvature'].values,alternative = 'two-sided')
			edge_pvals.append(p)
			edge_res['Edge'].append(e)
			edge_res['U'].append(u)
			edge_res['pval'].append(p)

		adjusted_pvals = mt.multipletests(edge_pvals,method = 'fdr_bh')
		edge_res['adj_pvals'] = adjusted_pvals[1]
		edge_results = pd.DataFrame(edge_res)
		edge_results = edge_results[edge_results['adj_pvals']<self.pvalue_thresh]
		edges = edge_results['Edge'].values
		signif_nodes = []
		for e in edges:
			genes = [x for x in e.split(";")]
			signif_nodes.extend(genes)
			
			
		gene_names = pd.unique(signif_nodes)
		col_indices = [self.gene_to_col[g] for g in gene_names]
		return gene_names, col_indices



	def compute_node_pvalues(self):
		node_data = pd.DataFrame(self.node_curvatures)
		node_res = defaultdict(list)
		node_pvals = []
		
		for n in pd.unique(node_data['Gene']):
			tempR = node_data[node_data['Gene'] == n]
			tempR = tempR[tempR['Response'] == 1]
			tempNR = node_data[node_data['Gene']==n]
			tempNR = tempNR[tempNR['Response']==0]
			u, p = mannwhitneyu(tempR['Raw'].values,tempNR['Raw'].values,alternative = 'two-sided')
			node_pvals.append(p)
			node_res['Gene'].append(n)
			node_res['U'].append(u)
			node_res['pval'].append(p)
		adjusted_pvals = mt.multipletests(node_pvals,method = 'fdr_bh')
		node_res['adj_pvals'] = adjusted_pvals[1]
		node_results = pd.DataFrame(node_res)
		node_results = node_results[node_results['adj_pvals']<self.pvalue_thresh]
		genes = pd.unique(node_results['Gene'].values)
		idxs = [self.gene_to_col[g] for g in genes]
		return genes, idxs

	def compute_normalized_node_pvalues(self):
		node_data = pd.DataFrame(self.node_curvatures)
		node_res = defaultdict(list)
		node_pvals = []
		
		for n in pd.unique(node_data['Gene']):
			tempR = node_data[node_data['Gene'] == n]
			tempR = tempR[tempR['Response'] == 1]
			tempNR = node_data[node_data['Gene']==n]
			tempNR = tempNR[tempNR['Response']==0]
			u, p = mannwhitneyu(tempR['Raw'].values,tempNR['Normalized'].values,alternative = 'two-sided')
			node_pvals.append(p)
			node_res['Gene'].append(n)
			node_res['U'].append(u)
			node_res['pval'].append(p)
		adjusted_pvals = mt.multipletests(node_pvals,method = 'fdr_bh')
		node_res['adj_pvals'] = adjusted_pvals[1]
		node_results = pd.DataFrame(node_res)
		node_results = node_results[node_results['adj_pvals']<self.pvalue_thresh]
		genes = pd.unique(node_results['Gene'].values)
		idxs = [self.gene_to_col[g] for g in genes]
		return genes, idxs
		
	