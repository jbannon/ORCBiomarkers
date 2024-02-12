from sklearn.metrics import confusion_matrix, roc_auc_score, auc, precision_recall_curve, accuracy_score
from sklearn.model_selection import train_test_split, GridSearchCV
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
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
import biomarkers as bm






def main(
	Config:Dict
	) -> None:
	
	
	
	
	drug, tissues, rng_seed, alpha, do_hops,num_hops, fdr_thresh,  min_TPM, num_iters, prob_thresh, train_pct, exp_type = \
		utils.unpack_parameters(config['EXPERIMENT_PARAMS']) 
	

	data_dir, geneset_dir, network_dir,result_dir = utils.unpack_parameters(config['DIRECTORIES'])
	
	path_method, min_degree_base, min_degree_exp, min_distance_base, min_distance_exp, sinkhorn_thresh, epsilon = \
		utils.unpack_parameters(config['CURVATURE_PARAMS'])

	pvalue_thresh, lcc_only = utils.unpack_parameters(config['FEATURE_PARAMS'])
	
	topology, weighting = utils.unpack_parameters(config['NETWORK_PARAMS'])
	

	weight_field, rename_field, density_field, edge_curvature_field, node_curve_field, norm_node_curve_field = \
		utils.unpack_parameters(config['FIELD_NAMES'])


	model_name, min_C, max_C, C_step, max_iter = utils.unpack_parameters(config["MODEL_PARAMS"])
	min_degree = min_degree_base**min_degree_exp
	min_distance = min_distance_base**min_distance_exp
	
	print(density_field)
	print(lcc_only)

	rng = np.random.RandomState(rng_seed)
	
	
	dataset_string = utils.DRUG_DATASET_MAP[drug]
	drug_targets = utils.fetch_drug_targets(drug)

	for tissue in tissues:
		expression_file = utils.make_file_path(data_dir,[dataset_string,drug, tissue],'expression','.csv')
		response_file = utils.make_file_path(data_dir,[dataset_string,drug,tissue],'response','.csv')
		diff_exp_file = utils.make_file_path(geneset_dir,[],"{d}_{t}_DE".format(d=drug, t=tissue),".csv")

		DE = pd.read_csv(diff_exp_file)
		DE = DE[DE['Thresh.Value']==fdr_thresh]
		
		response = pd.read_csv(response_file)		
		expression = pd.read_csv(expression_file)
		
		
		res_path = "".join([result_dir,"/".join(["classification",drug,"target_spread",tissue]),"/"])
		# res_path = "".join([result_dir,"/".join(["classification",exp_type,drug,tissue,lcc_string]),"/"])
		os.makedirs(res_path,exist_ok = True)
		
		
		res_name = "{p}monte_carlo.csv".format(p=res_path)
		stat_name = "{p}stats.txt".format(p=res_path)
		
		network_file = utils.make_file_path(network_dir,[dataset_string,topology],weighting,".pickle")
		with open(network_file,"rb") as istream:
			PPI_Graph = pickle.load(istream)
		
		common_genes = [x for x in drug_targets]
		for hop in range(1):
			# print("doing one hop")
			seeds = [x for x in common_genes if x in PPI_Graph.nodes()]
			
			one_hop = [x for x in seeds]
			
			for n in seeds:
				nbrs = PPI_Graph.neighbors(n)
				one_hop.extend([x for x in nbrs])
			common_genes = [x for x in one_hop]
			H = PPI_Graph.subgraph(common_genes)
			
		

		LCC_Graph = utils.harmonize_graph_and_geneset(PPI_Graph,common_genes,lcc_only)
		
		stat_string = "\n\tGraph has {n} nodes and {e} edges".format(n = len(LCC_Graph.nodes), e = len(LCC_Graph.edges))
		

		
		
		keep_cols = [x for x in LCC_Graph.nodes()]
		

		# subset_expression = expression[['Run_ID']+[x for x in LCC_Graph.nodes()]]
		subset_expression = expression[['Run_ID']+keep_cols]
		
		LCC_Graph, gene_to_idx, idx_to_gene = utils.rename_nodes(LCC_Graph,rename_field)
		subset_expression.set_index('Run_ID',inplace = True)

		
		gene_to_col = {}
		
		
		colnames = list(subset_expression.columns)
		
		for i in range(len(colnames)):
			gene_to_col[colnames[i]] = i
		
		X_ = subset_expression.values
		X_ = np.log2(X_+1)
		y = response['Response'].values
		
		
		model, param_grid = utils.make_model_and_param_grid(model_name,min_C,max_C,C_step,max_iter)
		
		results = defaultdict(list)

		feature_selector = bm.CurvatureFeatureSelector(
						base_graph = LCC_Graph,
						X = X_,
						y = y,
						gene_to_col = gene_to_col,
						gene_to_node_num = gene_to_idx,
						node_num_to_gene = idx_to_gene,
						alpha = alpha,
						pvalue_thresh = pvalue_thresh,
						weight_field = weight_field,
						curvature_field = edge_curvature_field,
						node_field = node_curve_field,
						norm_node_field = norm_node_curve_field,
						measure_name = density_field,
						min_distance = min_distance,
						min_degree = min_degree,
						sinkhorn_thresh = sinkhorn_thresh,
						epsilon = epsilon)

		feature_selector.compute_and_store_node_curvatures([gene_to_idx[x] for x in drug_targets])
		

		for feature in ['targets','node','norm-node']:
			
			if feature == 'targets':
				X = X_[:,[gene_to_col[x] for x in drug_targets]]
			elif feature == "node":
				X = feature_selector.X_node
			elif feature == "norm-node":
				X = feature_selector.X_norm
				
			
			for i in range(num_iters):
				X_train, X_test, y_train, y_test = train_test_split(X,y, train_size = train_pct, 
					random_state = rng,shuffle = True, stratify = y)
					
				clf = GridSearchCV(model, param_grid)
				clf.fit(X_train,y_train)
			
				train_preds_bin = clf.predict(X_train)
				train_preds_prob = clf.predict_proba(X_train)

				test_preds_bin = clf.predict(X_test)
				test_preds_prob = clf.predict_proba(X_test)
				train_acc = accuracy_score(y_train, train_preds_bin)
				test_acc = accuracy_score(y_test, test_preds_bin) 

				train_roc = roc_auc_score(y_train,train_preds_prob[:,1])
				test_roc = roc_auc_score(y_test,test_preds_prob[:,1])

				tn, fp, fn, tp = confusion_matrix(y_test, test_preds_bin,labels = [0,1]).ravel()
							

				results['drug'].append(drug)
				results['tissue'].append(tissue)
				results['feature'].append(feature)

				results['iter'].append(i)
				results['Train Accuracy'].append(train_acc)
				results['Train ROC_AUC'].append(train_roc)
				results['Test Accuracy'].append(test_acc)
				results['Test ROC_AUC'].append(test_roc)
				results['Test TN'].append(tn)
				results['Test FP'].append(fp)
				results['Test FN'].append(fn)
				results['Test TP'].append(tp)
	
		results = pd.DataFrame(results)
		
		
		results.to_csv(res_name)
		with open(stat_name,"w") as ostream:
			ostream.write(stat_string)




		
	
		
					


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)
	
	