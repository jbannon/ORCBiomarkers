import os
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
import itertools as it
from sklearn.model_selection import LeaveOneOut


def main(
	config:Dict
	) -> None:
	
	drug:str 
	rng_seed:int

	pval_thresh: float
	
	
	
	num_iters:int
	c_min:float
	c_max:float
	c_step:float
	model_types:List[str] 
	model_max_iters:int 

	drug, sources, targets, rng_seed, n_iters, exp_types = utils.unpack_parameters(config['EXPERIMENT_PARAMS'])
	data_dir, geneset_dir, network_dir,result_dir = utils.unpack_parameters(config['DIRECTORIES'])

	num_iters, c_min, c_max, c_step, model_max_iters = \
		utils.unpack_parameters(config['MODEL_PARAMS'])
	model_name, min_C, max_C, C_step, max_iter = utils.unpack_parameters(config["MODEL_PARAMS"])

	edge_pvs, node_pvs, norm_node_pvs = utils.unpack_parameters(config['FEATURE_PARAMS'])


	drug_targets = utils.fetch_drug_targets(drug)
	dataset_string = utils.DRUG_DATASET_MAP[drug]
	
	ctype_map = {'edge':'edge','node':'node','norm-node':'normalized_node'}
	pval_map = {'edge':edge_pvs,'node':node_pvs,'norm-node':norm_node_pvs}
	
	model, param_grid = utils.make_model_and_param_grid(model_name,min_C,max_C,C_step,max_iter)
	clf = GridSearchCV(model, param_grid)
	rng = np.random.RandomState(rng_seed)
	for source, target in tqdm.tqdm(it.product(sources,targets)):
		if source == target:
			continue
		res_path = "../results/transfer/{d}/".format(d=drug)
		os.makedirs(res_path, exist_ok = True)

		expression_file = utils.make_file_path(data_dir,[dataset_string,drug, target],'expression','.csv')
		response_file = utils.make_file_path(data_dir,[dataset_string,drug,target],'response','.csv')
		response = pd.read_csv(response_file)	
		expression = pd.read_csv(expression_file)
		results = defaultdict(list)
		for experiment_type in exp_types:
			print("working on: {t}".format(t=experiment_type))
			
			pval_dir = "../results/biomarkers/{et}/{d}/{s}/".format(et= experiment_type,d=drug, s = source)
			
			for ct in ['edge','node','norm-node']:
				pvalue_file = "{p}{c}_curvature_pvals.csv".format(p = pval_dir,c=ctype_map[ct])
				pvalues = pd.read_csv(pvalue_file,index_col = 0)
				for pv in pval_map[ct]:
					subset = pvalues[pvalues['adj_pvals']<pv]
					if ct in ['node','norm-node']:
						genes = list(pd.unique(subset['Gene']))
					else:
						genes = []
						for e in subset['Edge'].values:
							endpoints = [x for x in e.split(";")]
							genes.extend(endpoints)
						genes = [x for x in pd.unique(genes)]

					avoid_sampling = [x for x in genes]
					possible_samples = [i for i in expression.columns[1:] if i not in avoid_sampling]
					random_genes = np.random.choice(possible_samples,len(genes),replace = False)
					p = [x for x in possible_samples if x in genes]

					for feat_name, gene_names in zip([ct,ct+"-matched-random"],[genes,random_genes]):
						

						X = np.log2(expression[gene_names].values+1)
						y = response['Response'].values
						


						for i in tqdm.tqdm(range(n_iters)):
							X_train, X_test, y_train, y_test = train_test_split(X,y, test_size = 0.25, random_state = rng)
							
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
							results['source tissue'].append(source)
							results['target tissue'].append(target)
							results['type'].append(experiment_type)
							results['feature'].append(feat_name)
							results['pval'].append(pv)
							results['iter'].append(i)
							results['feature_dim'].append(len(gene_names))
							results['Train Accuracy'].append(train_acc)
							results['Train ROC_AUC'].append(train_roc)
							results['Test Accuracy'].append(test_acc)
							results['Test ROC_AUC'].append(test_roc)
							results['Test TN'].append(tn)
							results['Test FP'].append(fp)
							results['Test FN'].append(fn)
							results['Test TP'].append(tp)
		
		df = pd.DataFrame(results)
		df.to_csv("{p}{s}_to_{t}.csv".format(p=res_path, s=source, t=target),index = False)
			
			
		# gene_names = pd.unique(signif_nodes)
					

		










if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)