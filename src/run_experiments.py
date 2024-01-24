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
import yaml 
import utils
import pickle
import NetworkCurvature as nc
import statsmodels.stats.multitest as mt
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.linear_model import LogisticRegression
import tqdm
from sklearn.svm import LinearSVC
from sklearn.metrics import confusion_matrix, roc_auc_score, auc, precision_recall_curve, accuracy_score

def main(
	config:Dict
	) -> None:
	
	drug:str 
	rng_seed:int

	biomarker_types: List[str]
	genesets: List[str]
	topologies: List[str]
	cutoffs: List[Union[str,int]]

	
	edge_pval_thresh:float
	node_pval_thresh:float
	norm_node_pval_thresh:float
	
	
	
	num_iters:int
	c_min:float
	c_max:float
	c_step:float
	model_types:List[str] 
	model_max_iters:int 
	
	drug, rng_seed, biomarker_types, genesets, topologies, cutoffs, edge_pval_thresh, node_pval_thresh, norm_node_pval_thresh = \
		utils.unpack_parameters(config['EXPERIMENT_PARAMS'])


	num_iters, c_min, c_max, c_step, model_types, model_max_iters = \
		utils.unpack_parameters(config['CROSSVAL_PARAMS'])

	tissues = utils.DRUG_TISSUE_MAP[drug]
	np.random.seed(rng_seed)

	
	ds = utils.DRUG_DATASET_MAP[drug]
	targets = utils.fetch_drug_targets(drug)
	

	for tissue in tissues:
		expr_name = "../data/expression/{ds}/{d}/{t}/expression.csv".format(ds=ds,d=drug,t=tissue)
		resp_name = "../data/expression/{ds}/{d}/{t}/response.csv".format(ds= ds,d=drug,t=tissue)
		expr = pd.read_csv(expr_name, index_col = 0)
		resp = pd.read_csv(resp_name, index_col = 0)
		results = defaultdict(list)
		

		res_name = "../results/experiments/{d}_{t}_results.csv".format(d=drug, t = tissue)



		for b_type in biomarker_types:
			for geneset in genesets:
				for topology in topologies:
					for cutoff in cutoffs:
						if drug == 'Atezo' and tissue == 'KIRC' and  geneset == 'DE':
							continue
						edge_fname = utils.make_file_path(
							base_dir = "../results/biomarkers/",
							path_names = [b_type, geneset, drug, tissue, str(cutoff), topology, "NoHop"],
							fname = "edge_curvature_pvals",
							ext = ".csv")

						node_fname = utils.make_file_path(
							base_dir = "../results/biomarkers/",
							path_names = [b_type, geneset, drug, tissue, str(cutoff), topology, "NoHop"],
							fname = "node_curvature_pvals",
							ext = ".csv")

						norm_node_fname = utils.make_file_path(
							base_dir = "../results/biomarkers/",
							path_names = [b_type, geneset, drug, tissue, str(cutoff), topology, "NoHop"],
							fname = "normalized_node_curvature_pvals",
							ext = ".csv")
						


						edge_pvs = pd.read_csv(edge_fname,index_col=0)
						node_pvs = pd.read_csv(node_fname,index_col=0)
						norm_node_pvs = pd.read_csv(norm_node_fname,index_col=0)


						edge_genes = utils.process_pvalue_data(
							pval_df = edge_pvs,
							gene_col = "Edge",
							pval_thresh = edge_pval_thresh,
							pval_col = 'adj_pvals')


						num_edge_genes = len(edge_genes)

						
						node_genes = utils.process_pvalue_data(
							pval_df = node_pvs,
							gene_col = "Gene",
							pval_thresh = node_pval_thresh,
							pval_col = 'adj_pvals')

						num_node_genes = len(node_genes)
						

						norm_node_genes =  utils.process_pvalue_data(
							pval_df = norm_node_pvs,
							gene_col = "Gene",
							pval_thresh = norm_node_pval_thresh,
							pval_col = 'adj_pvals')

						num_norm_node_genes = len(norm_node_genes)


						edge_randoms = np.random.choice(expr.columns[1:],num_edge_genes,replace=False)
						node_randoms = np.random.choice(expr.columns[1:],num_node_genes,replace=False)
						norm_node_randoms = np.random.choice(expr.columns[1:],num_norm_node_genes,replace=False)


						
						for name, feature_genes in zip(
							["targets","edge","node", "norm-node","edge-random","node-random","norm-node-random"],
							[targets, edge_genes, node_genes, norm_node_genes, edge_randoms, node_randoms, norm_node_randoms]):
							
							
							X = np.log2(expr[feature_genes].values+1)
							y = resp['Response'].values
							
							for model_name in tqdm.tqdm(model_types):
								model, param_grid = utils.make_model_and_param_grid(
									model_name = model_name,
									reg_min = c_min,
									reg_max = c_max,
									reg_step = c_step,
									model_max_iters = model_max_iters
									)
								

								
								for k in tqdm.tqdm(range(num_iters),leave=False):
									X_train, X_test, y_train, y_test = train_test_split(X,y,stratify = y)
									

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
									results['biomarker'].append(b_type)
									results['geneset'].append(geneset)
									results['topology'].append(topology)
									results['cutoff'].append(cutoff)
									results['model'].append(model_name)

									results['iter'].append(k)
									results['model'].append(model_name)
									results['geneset'].append(name)
									results['num_features'].append(len(geneset))
									results['Train Accuracy'].append(train_acc)
									results['Train ROC_AUC'].append(train_roc)
									results['Test Accuracy'].append(test_acc)
									results['Test ROC_AUC'].append(test_roc)
									results['Test TN'].append(tn)
									results['Test FP'].append(fp)
									results['Test FN'].append(fn)
									results['Test TP'].append(tp)

		df = pd.DataFrame(results)
		df.to_csv(res_name,index = False)

						








if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)