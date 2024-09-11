from sklearn.feature_selection import SelectKBest, chi2, f_classif, mutual_info_classif
import warnings
warnings.filterwarnings("ignore")
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



def select_biomarker_subset(
	genes:List[str],
	expression:pd.DataFrame,
	response:pd.DataFrame,
	seed_genes:List[str],
	thresh:float = 0.1
	):
	
	expr_ = expression[['Run_ID']+genes]	
	for g in genes:
		expr_[g]=np.log2(expr_[g]+1)

		

	data = expr_.merge(response, on='Run_ID')
		
			
	data.drop(columns = ['Run_ID'],inplace=True)
	data = data.melt(id_vars=['Response'],var_name = 'Gene')
	
	p_vals = []
	adj_pvals = []
	_g = []
	# print(data)
	for gene in pd.unique(data['Gene']):
		temp = data[data['Gene']==gene]
		temp_r =temp[temp['Response']==1]
		temp_nr =temp[temp['Response']==0]
		r_exp = temp_r['value'].values
		nr_exp = temp_nr['value'].values
		U,p = mannwhitneyu(r_exp, nr_exp)
		p_vals.append(p)
		_g.append(gene)

	adj_pvals = mt.multipletests(p_vals,method = 'fdr_bh')

	_kept_genes = [_g[i] for i in np.where(adj_pvals[1]<=thresh)[0]]
	_kept_genes = [x for x in _kept_genes if x not in seed_genes]
	
	return _kept_genes

def main(
	config:Dict
	) -> None:
	
	source_drug:str 
	target_drug:str
	rng_seed:int

	pval_thresh: float
	
	
	# num_perms:int = 100
	num_iters:int
	c_min:float
	c_max:float
	c_step:float
	model_types:List[str] 
	model_max_iters:int 

	source_drug, target_drug, source_tissues, target_tissues, rng_seed, n_iters = \
		utils.unpack_parameters(config['EXPERIMENT_PARAMS'])
	

	data_dir, geneset_dir, network_dir,result_dir = \
		utils.unpack_parameters(config['DIRECTORIES'])

	num_iters, c_min, c_max, c_step, model_max_iters = \
		utils.unpack_parameters(config['MODEL_PARAMS'])
	
	model_name, min_C, max_C, C_step, max_iter = \
		utils.unpack_parameters(config["MODEL_PARAMS"])

	edge_pvs, node_pvs, norm_node_pvs = \
		utils.unpack_parameters(config['FEATURE_PARAMS'])


	drug_targets = utils.fetch_drug_targets(target_drug)
	target_dataset_string = utils.DRUG_DATASET_MAP[target_drug]
	source_dataset_string = utils.DRUG_DATASET_MAP[source_drug]
	
	ctype_map = {'edge':'edge','node':'node','norm-node':'normalized_node'}
	pval_map = {'edge':edge_pvs,'node':node_pvs,'norm-node':norm_node_pvs}
	
	model, param_grid = utils.make_model_and_param_grid(model_name,min_C,max_C,C_step,max_iter)
	clf = GridSearchCV(model, param_grid)
	

	rng = np.random.RandomState(rng_seed)
	
	with open("../data/genesets/kegg.txt", "r") as istream:
		lines = istream.readlines()

	kegg_genes = [x.rstrip() for x in lines]

	
	gene_lens = []
	for source_tissue, target_tissue in it.product(source_tissues,target_tissues):
		if source_tissue == target_tissue and source_drug == target_drug:
			continue

		Filtered_DE_file = f"../results/biomarkers/{source_drug}/{source_tissue}/Filtered_DE_genes.txt"

		with open(Filtered_DE_file, "r") as istream:
			Filtered_DE_genes = istream.readlines()
			Filtered_DE_genes = [x.rstrip() for x in Filtered_DE_genes]
		

		expression_file = utils.make_file_path(
			data_dir,
			[target_dataset_string,target_drug,
			target_tissue],
			'expression','.csv'
			)
		response_file = utils.make_file_path(
			data_dir,[
			target_dataset_string,target_drug,
			target_tissue],
			'response','.csv')

		target_response = pd.read_csv(response_file)	
		target_expression = pd.read_csv(expression_file)

		expression_file = utils.make_file_path(
			data_dir,
			[source_dataset_string,source_drug,source_tissue],
			'expression','.csv')

		response_file = utils.make_file_path(
			data_dir,
			[source_dataset_string,source_drug,source_tissue],
			'response','.csv')
		
		
		source_response = pd.read_csv(response_file)	
		source_expression = pd.read_csv(expression_file)



		# _kegg = [x for x in kegg_genes if x in expression.columns]
		# gene_lens.append(('kegg',len(_kegg)))

		# res_path = "../results/transfer/{d}/monte_carlo/".format(d=drug)
		res_path = f"../results/transfer/{source_drug}/{source_tissue}/{target_drug}/{target_tissue}/"
		os.makedirs(res_path, exist_ok = True)
		results = defaultdict(list)
		
			
		pval_dir = f"../results/biomarkers/{source_drug}/{source_tissue}/"
			
		for ct in ['edge','node','norm-node']:
		# for ct in ['edge']:	
			pvalue_file = f"{pval_dir}{ctype_map[ct]}_curvature_pvals.csv"
			pvalues = pd.read_csv(pvalue_file,index_col = 0)
			


			for pv in tqdm.tqdm(pval_map[ct]):

				results = defaultdict(list)
				subset = pvalues[pvalues['adj_pvals']<pv]
				if subset.shape[0]==0:
					continue
				if ct in ['node','norm-node']:
					genes = list(pd.unique(subset['Gene']))
				else:
					genes = []
					for e in subset['Edge'].values:
						endpoints = [x for x in e.split(";")]
						genes.extend(endpoints)
					genes = [x for x in pd.unique(genes)]
				
				
				filtered_genes = select_biomarker_subset(
					genes,
					source_expression,
					source_response,
					Filtered_DE_genes)


				# print(filtered_genes)
				
				# avoid_sampling = [x for x in filtered_genes]

				# possible_samples = [i for i in target_expression.columns[1:] if i not in avoid_sampling]
				
				# random_genes = np.random.choice(possible_samples,len(genes),replace = False)
				results = []
				for gs_name, geneset in zip([
					'filtered','all'],[filtered_genes,genes]):
					results.append(gs_name.upper()+"\n--------\n")
					
					for g in geneset:
						target_expression[g]=np.log2(target_expression[g]+1)
				
					X = target_expression[filtered_genes]
					# np.log2(target_expression[[x for x in filtered_genes if x in target_expression.columns]].values+1)
					y = target_response['Response'].values
					# int(np.floor(X.shape[1]/2))
				


					k = min(int(np.floor(X.shape[1]/2)),10)
					all_names =[]
					counting_names = []
					for name, meas in zip(['chi2','f_meas','mi'],
						[chi2, f_classif, mutual_info_classif]):
						selector = SelectKBest(meas, k = k)
						X_new = selector.fit_transform(X, y)
						gene_names = selector.get_feature_names_out()
						all_names.append(gene_names)
						counting_names.extend(gene_names)
						results.append(f"For measure {name} the {k} best are: {gene_names}\n")

					common_genes = set.intersection(*[set(x) for x in all_names])
					results.append(f"\nThe intersection of all gives: {common_genes}\n")
					all_names = [x for name_list in all_names for x in name_list]
					counts = pd.value_counts(counting_names)
					two = list(counts.iloc[np.where(counts>=2)[0],].index)
					results.append(f"The two count: {two}\n")
					results.append("\n\n")


				with open(f"{res_path}{ct}_{pv}.txt","w") as ostream:
					ostream.writelines(results)

		










if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)