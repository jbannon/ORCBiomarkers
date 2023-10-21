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
import statsmodels.stats.multitest.multipletests






def main(
	Config:Dict
	) -> None:
	
	# drug:str
	# tissue:str
	# alphas:List[float] 
	# num_iters:int
	# expr_dir:str
	# geneset_dir:str
	# do_one_hop:bool
	# norm: str
	# rng_seed: int

	# fdr_threshold: float
	# count_cutoff: int


	np.random.seed(12345)
	# drug, tissue, alphas, num_iters, expr_dir, geneset_dir, do_one_hop, norm =\
	# 	utils.UNPACK_PARAMETERS(config['EXPERIMENT_PARMS']) 
	
	# fdr_threshold, count_cutoff = utils.UNPACK_PARAMETERS(config['DE_PARAMS'])



	utils.make_file_name(['a','b','c'])
	expr_file = "{data_directory}/{drug}/{tissue}/expression.csv".format()
	expr = pd.read_csv('../data/expression/Atezo/KIRC/expression.csv')
	resp = pd.read_csv('../data/expression/Atezo/KIRC/response.csv')
	DE = pd.read_csv("../data/genesets/Atezo_KIRC_DE.csv",index_col = 0)
	



	DE = DE[DE['Thresh.Value']==fdr_threshold]
	DE = DE[DE['Count']>=count_cutoff]
	
	DE_genes = list(DE['Gene'].values)
	
	if do_one_hop:
		pass
		# collect the neighbors then reduce to the LCC

	responder_IDs = resp[resp['Response']==1]['Run_ID'].tolist()
	nonResponder_IDs = resp[resp['Response']==0]['Run_ID'].tolist()

	num_samples = np.round(min(len(responder_IDs)/2,len(nonResponder_IDs)/2))

	for i in range(num_iters):
		pass 
		# sample_responders
		# sample non
		responder_samples = np.random.choice(responder_IDs, num_samples, replace = False)
		nonResponder_samples = np.random.choice(nonResponder_IDs, num_samples, replace = False)

	expr = expr[['Run_ID']+DE_genes]
	

	responder_measurements = expr[expr['Run_ID'].isin(responder_samples)]
	nonResponder_measurements = expr[expr['Run_ID'].isin(nonResponder_samples)]

	G_0 = build_network(nonResponder_measurements)
	print(nx.is_connected(G_0))
	assign_densities(G_0)



if __name__ == '__main__':
	# parser = argparse.ArgumentParser()
	# parser.add_argument("-config",help="The config file for these experiments")
	# args = parser.parse_args()
	
	# with open(args.config) as file:
	# 	config = yaml.safe_load(file)

	# main(config)
	main({})
	