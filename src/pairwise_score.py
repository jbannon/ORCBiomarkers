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

	drug, sources, targets, rng_seed, n_iters = utils.unpack_parameters(config['EXPERIMENT_PARAMS'])
	data_dir, geneset_dir, network_dir,result_dir = utils.unpack_parameters(config['DIRECTORIES'])

	num_iters, c_min, c_max, c_step, model_max_iters = \
		utils.unpack_parameters(config['MODEL_PARAMS'])
	model_name, min_C, max_C, C_step, max_iter = utils.unpack_parameters(config["MODEL_PARAMS"])

	drug_targets = utils.fetch_drug_targets(drug)

	dataset_string = utils.DRUG_DATASET_MAP[drug]

	
	for tissue in targets:
		edge_file = "../results/biomarkers/{d}/{t}/edge_curvature_pvals.csv".format(d=drug,t=tissue)
		edges = pd.read_csv(edge_file,index_col = 0)
		expression_file = utils.make_file_path(data_dir,[dataset_string,drug, tissue],'expression','.csv')
		response_file = utils.make_file_path(data_dir,[dataset_string,drug,tissue],'response','.csv')
		response = pd.read_csv(response_file)	
		expression_file = pd.read_csv(expression_file)
		edges = edges[edges['adj_pvals']<=0.005]
		edges = edges['Edge'].values

		for e in edges:
			print(e)




if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)
