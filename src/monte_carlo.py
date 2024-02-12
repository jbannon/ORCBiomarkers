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
	

	rng = np.random.RandomState(rng_seed)
	
	
	dataset_string = utils.DRUG_DATASET_MAP[drug]
	drug_targets = utils.fetch_drug_targets(drug)
	lcc_string = "lcc_only" if lcc_only else "full_graph"

	for tissue in tissues:
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
		
		
		
		DE_gene_list = utils.empirical_bayes_gene_selection(DE,prob_thresh)
		stat_string = "\n\t{d} - {t}\n".format(d=drug, t= tissue)
		stat_string += "\n\t{k} DE genes after filtering".format(k=len(DE_gene_list))
		stat_string += "\n\t{p} patients in the dataset".format(p = response.shape[0])

		# genes that are differentially expressed but also greater than tpm threshold in all entries
		common_genes = [x for x in DE_gene_list if x in keep_genes]
		DE_genes = [x for x in common_genes]

		res_path = "".join([result_dir,"/".join(["classification",drug,lcc_string,tissue]),"/"])
		#res_path = "".join([result_dir,"/".join(["classification",exp_type,drug,tissue,lcc_string]),"/"])
		os.makedirs(res_path,exist_ok = True)
		
		lcc_string = "_lcc" if lcc_only else ""
		res_name = "{p}monte_carlo.csv".format(p=res_path,l = lcc_string)
		stat_name = "{p}stats.txt".format(p=res_path, l=lcc_string)
		
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
		
		# should be fine if we do a pandas unique for the LCC nodes, drug targets, and pre-common

		LCC_Graph = utils.harmonize_graph_and_geneset(PPI_Graph,common_genes,lcc_only)
		
		stat_string += "\n\tGraph has {n} nodes and {e} edges".format(n = len(LCC_Graph.nodes), e = len(LCC_Graph.edges))
		

		
		keep_cols = [x for x in list(pd.unique([g for g in LCC_Graph.nodes()]+DE_genes+drug_targets))]
		

		# subset_expression = expression[['Run_ID']+[x for x in LCC_Graph.nodes()]]
		subset_expression = expression[['Run_ID']+keep_cols]
		
		LCC_Graph, gene_to_idx, idx_to_gene = utils.rename_nodes(LCC_Graph,rename_field)
		subset_expression.set_index('Run_ID',inplace = True)

		
		gene_to_col = {}
		
		
		colnames = list(subset_expression.columns)
		
		for i in range(len(colnames)):
			gene_to_col[colnames[i]] = i
		
		X = subset_expression.values
		X = np.log2(X+1)
		y = response['Response'].values
		
		
		model, param_grid = utils.make_model_and_param_grid(model_name,min_C,max_C,C_step,max_iter)
		
		results = defaultdict(list)
		
	
		
		for i in range(num_iters):
			X_train, X_test, y_train, y_test = train_test_split(X,y, train_size = train_pct, 
				random_state = rng,shuffle = True, stratify = y)

			feature_selector = bm.CurvatureFeatureSelector(
						base_graph = LCC_Graph,
						X = X_train,
						y = y_train,
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


			feature_selector.compute_curvatures()

			for feature in ['targets','DE','edge','node',
				'node-norm',
				'edge+targets','node+targets',
				'node-norm+targets']:
				

				
				if feature == 'targets':
					names = ['targets']
					indices = [[gene_to_col[x] for x in drug_targets]]
				elif feature == 'DE':
					names = ["DiffExp"]
					indices = [[gene_to_col[x] for x in DE_genes]]
				else:
					if feature[-8:]=="+targets":
						feat_type = feature[:-8]
					else:
						feat_type = feature
					
					genes, idxs  = feature_selector.get_features(feature_type = feat_type)
					if feature[-8:]=="+targets":
						idxs = idxs + [gene_to_col[x] for x in drug_targets]

					avoid_sampling = [x for x in idxs] + [gene_to_col[x] for x in drug_targets]
					possible_samples = [i for i in np.arange(X.shape[1]) if i not in avoid_sampling]
					
					random_genes = np.random.choice(np.arange(X.shape[1]),len(idxs),replace = False)
					names = [feature,feature+"-random"]
					indices = [idxs,random_genes]
				
				for feature_name, col_idxs in zip(names, indices):
					if len(col_idxs)==0:
						results['drug'].append(drug)
						results['tissue'].append(tissue)
						results['feature'].append(feature_name)

						results['iter'].append(i)
						results['feature_dim'].append(-1)
						results['Train Accuracy'].append(-1)
						results['Train ROC_AUC'].append(-1)
						results['Test Accuracy'].append(-1)
						results['Test ROC_AUC'].append(-1)
						results['Test TN'].append(-1)
						results['Test FP'].append(-1)
						results['Test FN'].append(-1)
						results['Test TP'].append(-1)
		
					else:
						clf = GridSearchCV(model, param_grid)
						clf.fit(X_train[:,col_idxs],y_train)
					
						train_preds_bin = clf.predict(X_train[:,col_idxs])
						train_preds_prob = clf.predict_proba(X_train[:,col_idxs])

						test_preds_bin = clf.predict(X_test[:,col_idxs])
						test_preds_prob = clf.predict_proba(X_test[:,col_idxs])
						train_acc = accuracy_score(y_train, train_preds_bin)
						test_acc = accuracy_score(y_test, test_preds_bin) 

						train_roc = roc_auc_score(y_train,train_preds_prob[:,1])
						test_roc = roc_auc_score(y_test,test_preds_prob[:,1])

						tn, fp, fn, tp = confusion_matrix(y_test, test_preds_bin,labels = [0,1]).ravel()
									

						results['drug'].append(drug)
						results['tissue'].append(tissue)
						results['feature'].append(feature_name)

						results['iter'].append(i)
						results['feature_dim'].append(len(col_idxs))
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
	
	