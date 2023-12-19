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
# import statsmodels.stats.multitest.multipletests
import pickle
import NetworkCurvature as nc
import statsmodels.stats.multitest as mt
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.linear_model import LogisticRegression
import tqdm



def main(
	config:Dict
	) -> None:
	
	drug:str = 'Atezo'
	tissue:str = 'BLCA'
	# alphas:List[float] 
	alpha:float = 0.5
	num_iters:int
	expr_dir:str
	geneset_dir:str
	do_one_hop:bool
	norm: str
	rng_seed: int 
	fdr_threshold: float
	count_cutoff: int
	
	rng_seed:int = 12345
	np.random.seed(rng_seed)



	expr = pd.read_csv("../data/expression/cri/{d}/{t}/expression.csv".format(d=drug,t=tissue))
	resp = pd.read_csv("../data/expression/cri/{d}/{t}/response.csv".format(d=drug,t=tissue))
	

	edge_pvs = pd.read_csv("../results/biomarkers/correlation/COSMIC/Atezo/BLCA/EB/dense/NoHop/edge_curvature_pvals.csv")
	node_pvs = pd.read_csv("../results/biomarkers/correlation/COSMIC/Atezo/BLCA/EB/dense/NoHop/node_curvature_pvals.csv")
	norm_node_pvs = pd.read_csv("../results/biomarkers/correlation/COSMIC/Atezo/BLCA/EB/dense/NoHop/node_curvature_pvals.csv")



	thresh = 0.1

	edge_genes = edge_pvs[edge_pvs['adj_pvals']<=thresh]
	edges = edge_genes['Edge'].values
	edgeNodes = []
	for edge in edges:
		nodes = [x for x in edge.split(";")]
		for node in nodes:
			edgeNodes.append(node)

	edgeGenes = list(pd.unique(edgeNodes))
	results = defaultdict(list)
	random_genes = np.random.choice(expr.columns[1:],len(edgeGenes),replace=False)


	for balance in tqdm.tqdm([True,False]):
		for name, geneset in zip(["target","orc-edge",'random'],[['CD274'], 
			edgeGenes,
			random_genes]):
			
			X = np.log2(expr[geneset].values+1)
			y = resp['Response'].values

			param_grid = {'clf__C':np.arange(0.1,1,0.1)}
			
			if balance:
				cw = 'balanced'
			else:
				cw = None
			model = Pipeline([('preproc',StandardScaler()),('clf',LogisticRegression(class_weight = cw ))])

			for k in tqdm.tqdm(range(50),leave=False):
				X_train, X_test, y_train, y_test = train_test_split(X,y)
				

				clf = GridSearchCV(model, param_grid)
				clf.fit(X_train,y_train)


				preds = clf.best_estimator_.predict(X_test)
				
				results['iter'].append(k)
				results['geneset'].append(name)
				results['balanced'].append(str(balance))

				results['acc'].append(accuracy_score(y_test,preds))



	df = pd.DataFrame(results)
	for b in pd.unique(df['balanced']):
		temp = df[df['balanced']==b]
		sns.boxplot(x='geneset', y = 'acc',data = temp)
		plt.title("{d}\t{t}\t{bal}".format(d=drug,t=tissue,bal=b))
		plt.savefig("atezo_{b}.png".format(b=b))
		plt.close()








if __name__ == '__main__':
	main({})