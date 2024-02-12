import time 
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
import pickle
import NetworkCurvature as nc
import statsmodels.stats.multitest as mt
from scipy.stats import mannwhitneyu, wilcoxon
import tqdm



def main(config:Dict):
		
	drug, tissues, rng_seed, alpha, do_hops,num_hops, fdr_thresh,  min_TPM, num_iters, prob_thresh, train_pct, exp_type = \
		utils.unpack_parameters(config['EXPERIMENT_PARAMS']) 
	
	
	data_dir, geneset_dir, network_dir,result_dir = utils.unpack_parameters(config['DIRECTORIES'])
	
	path_method, min_degree_base, min_degree_exp, min_distance_base, min_distance_exp, sinkhorn_thresh, epsilon = \
		utils.unpack_parameters(config['CURVATURE_PARAMS'])

	pvalue_thresh, lcc_only = utils.unpack_parameters(config['FEATURE_PARAMS'])


	min_degree = min_degree_base**min_degree_exp
	min_distance = min_distance_base**min_distance_exp
	topology, weighting = utils.unpack_parameters(config['NETWORK_PARAMS'])
	
	weight_field, rename_field, density_field, edge_curvature_field, node_curve_field, norm_node_curve_field = \
		utils.unpack_parameters(config['FIELD_NAMES'])

	rng = np.random.default_rng(rng_seed)
	
	dataset_string = utils.DRUG_DATASET_MAP[drug]
	# drug_target = utils.TARGET_GENE_MAP[utils.DRUG_TARGET_MAP[drug]]
	lcc_string = "lcc_only" if lcc_only else "full_graph"
	for tissue in tissues:
		print('starting {t}'.format(t=tissue))
		expression_file = utils.make_file_path(data_dir,[dataset_string,drug, tissue],'expression','.csv')
		response_file = utils.make_file_path(data_dir,[dataset_string,drug,tissue],'response','.csv')
		diff_exp_file = utils.make_file_path(geneset_dir,[],"{d}_{t}_DE".format(d=drug, t=tissue),".csv")

		DE = pd.read_csv(diff_exp_file)
		DE = DE[DE['Thresh.Value']==fdr_thresh]

		response = pd.read_csv(response_file)	

		expression = pd.read_csv(expression_file)
		
		trimmed_expression = expression[expression.columns[1:]]
		
		
		trimmed_expression = np.sum(trimmed_expression>min_TPM,axis=0)/expression.shape[0]
		trimmed_expression = trimmed_expression[trimmed_expression>0.9]
		
		# trimmed_expression = (trimmed_expression>min_TPM).all(axis=0)
		# keep_genes = list(trimmed_expression[trimmed_expression].index)
		
		keep_genes = list(trimmed_expression.index)
		
		keep_cols = ['Run_ID'] + keep_genes
		gene_list = utils.empirical_bayes_gene_selection(DE,prob_thresh)
		
		expression = expression[keep_cols]
		
		
		
		DE_genes = [x for x in gene_list]
		num_DE = len(DE_genes)
		common_genes = [x for x in gene_list if x in expression.columns[1:]]
		filtered_DE = [x for x in common_genes]
		num_filtered_DE = len(filtered_DE)
		
		res_path = "".join([result_dir,"/".join(["biomarkers",drug,lcc_string,tissue]),"/"])
		#res_path = "".join([result_dir,"/".join(["biomarkers",exp_type,drug,tissue,lcc_string]),"/"])
		os.makedirs(res_path,exist_ok = True)

		with open(res_path+"DE_genes.txt","w") as ostream:
			ostream.writelines([x+"\n" for x in DE_genes])

		with open(res_path+"Filtered_DE_genes.txt","w") as ostream:
			ostream.writelines([x+"\n" for x in filtered_DE])
		

					
		network_file = utils.make_file_path(network_dir,[dataset_string,topology],weighting,".pickle")
		with open(network_file,"rb") as istream:
			PPI_Graph = pickle.load(istream)

		if do_hops:
			for hop in range(num_hops):
			# print("doing one hop")
				seeds = [x for x in common_genes if x in PPI_Graph.nodes()]
				one_hop = [x for x in seeds]
				for n in seeds:
					nbrs = PPI_Graph.neighbors(n)
					one_hop.extend([x for x in nbrs if x in keep_genes])
				common_genes = [x for x in one_hop]
		
		
		LCC_Graph = utils.harmonize_graph_and_geneset(PPI_Graph,common_genes,lcc_only)
		print(LCC_Graph)
		if len(LCC_Graph.edges)==0:
			with open(res_path+"empty_graph.txt","w") as ostream:
				ostream.writelines(["graph had no edges"])
			continue

		subset_expression = expression[['Run_ID']+[x for x in LCC_Graph.nodes()]]
		data = subset_expression.merge(response, on = 'Run_ID')	
		LCC_Graph, gene_to_idx, idx_to_gene = utils.rename_nodes(LCC_Graph,rename_field)


		edge_curvatures = defaultdict(list)
		node_curvatures = defaultdict(list)

		for index, row in tqdm.tqdm(data.iterrows(),total=data.shape[0],leave=False):
			G = LCC_Graph.copy()
			sample_response = row["Response"]
			
			for edge in G.edges():
				node1, node2 = edge[0], edge[1]
				gene1, gene2 = idx_to_gene[node1], idx_to_gene[node2]
				weight = np.round(np.log2(row[gene1]+1)*np.log2(row[gene2]+1),5)
				G[node1][node2]['weight'] = weight

			orc = nc.OllivierRicciCurvature(
				G= G,
				alpha = alpha, 
				weight_field = weight_field,
				path_method = path_method,
				curvature_field = edge_curvature_field,
				node_field = node_curve_field,
				norm_node_field = norm_node_curve_field,
				measure_name = density_field,
				min_distance = min_distance,
				min_degree = min_degree,
				sinkhorn_thresh = sinkhorn_thresh,
				epsilon = epsilon
				)
			
			orc.compute_curvatures()
			G = orc.G.copy()
		
		

			for edge in G.edges(data=True):
				edge_curvatures['Patient'].append(row['Run_ID'])
				edge_curvatures['Gene1'].append(idx_to_gene[edge[0]])
				edge_curvatures['Gene2'].append(idx_to_gene[edge[1]])
				edge_string= "{a};{b}".format(a=idx_to_gene[edge[0]] ,b = idx_to_gene[edge[1]])				
				edge_curvatures['Edge'].append(edge_string)
				edge_curvatures['Curvature'].append(edge[2][edge_curvature_field])
				edge_curvatures['Response'].append(sample_response)
			
			for node in G.nodes(data= True):
				# print(node)
				node_curvatures['patient'].append(row['Run_ID'])
				node_curvatures['Gene'].append(node[1][rename_field])
				node_curvatures['Raw'].append(node[1][node_curve_field])
				node_curvatures['Normalized'].append(node[1][norm_node_curve_field])
				node_curvatures['Response'].append(sample_response)
		
		edge_data = pd.DataFrame(edge_curvatures)
		node_data = pd.DataFrame(node_curvatures)
		
		edge_res = defaultdict(list)
		node_res = defaultdict(list)
		norm_res = defaultdict(list)
	
		node_pvals = []
		norm_pvals = []
		edge_pvals = []
		print("processing edge info")
		for e in tqdm.tqdm(pd.unique(edge_data['Edge'])):
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
		

	
		print("processing node info ")
		for n in tqdm.tqdm(pd.unique(node_data['Gene'])):
			tempR = node_data[node_data['Gene'] == n]
			tempR = tempR[tempR['Response'] == 1]
			tempNR = node_data[node_data['Gene']==n]
			tempNR = tempNR[tempNR['Response']==0]
			u, p = mannwhitneyu(tempR['Normalized'].values,tempNR['Normalized'].values,alternative = 'two-sided')
			
			norm_pvals.append(p)
			norm_res['Gene'].append(n)
			norm_res['U'].append(u)
			norm_res['pval'].append(p)

			u, p = mannwhitneyu(tempR['Raw'].values,tempNR['Raw'].values,alternative = 'two-sided')
			node_pvals.append(p)
			node_res['Gene'].append(n)
			node_res['U'].append(u)
			node_res['pval'].append(p)


		node_adjusted_pvals = mt.multipletests(node_pvals,method = 'fdr_bh')
		node_res['adj_pvals'] = node_adjusted_pvals[1]

		norm_node_adjusted_pvals = mt.multipletests(norm_pvals,method = 'fdr_bh')
		norm_res['adj_pvals'] = norm_node_adjusted_pvals[1]

		node_results = pd.DataFrame(node_res)
		norm_results = pd.DataFrame(norm_res)
		edge_results = pd.DataFrame(edge_res)
		norm_results.to_csv(res_path+"normalized_node_curvature_pvals.csv")
				
		node_results.to_csv(res_path+"node_curvature_pvals.csv")
				
		edge_results.to_csv(res_path+"edge_curvature_pvals.csv")
				
		edge_data.to_csv(res_path+"edge_curvatures.csv")
				
		node_data.to_csv(res_path+"node_curvatures.csv")
		
		stat_string = "alpha:\t{a}\nfdrThresh:\t{t}\nProb Thresh:\t{g}".format(a=alpha,t=fdr_thresh,g = prob_thresh)
		stat_string = stat_string + "\nMinimum TPM:\t{m}".format(m=min_TPM)
		stat_string = stat_string + "\nTake Hops?\t{m}".format(m=do_hops)
		stat_string = stat_string + "\nNum Hops:\t{s}".format(s=num_hops)
		stat_string = stat_string + "\nAlpha\t{e}".format(e=alpha)
		stat_string = stat_string + "\nLCC has {n} nodes and {e} edges".format(n=len(LCC_Graph.nodes),e=len(LCC_Graph.edges))
		stat_string = stat_string + "\nNum DE Genes Before Filtering:\t{n}".format(n=num_DE)
		stat_string = stat_string + "\nNum DE Genes After Filtering:\t{n}".format(n=num_filtered_DE)
		statname = res_path + "stats.txt"


		with open(statname,"w") as ostream:
			ostream.write(stat_string)

		

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)
	


	
				
				