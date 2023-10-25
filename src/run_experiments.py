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
	# connectivity: s/tr = 'sparse' 
	weighting: str = 'unweighted'
	np.random.seed(rng_seed)



	expr = pd.read_csv("../data/expression/{d}/{t}/expression.csv".format(d=drug,t=tissue))
	resp = pd.read_csv("../data/expression/{d}/{t}/response.csv".format(d=drug,t=tissue))
	DE = pd.read_csv("../data/genesets/{d}_{t}_DE.csv".format(d=drug,t=tissue),index_col = 0)



	DRUG_TARGET_MAP = {'Atezo':'PD-L1','Pembro':'PD1','Nivo':'PD1','Ipi':'CTLA4'}


	TARGET_GENE_MAP = {'PD-L1':'CD274', 'PD1':'PDCD1', 'CTLA4':'CTLA4'}


	results = defaultdict(list)
	random_genes = np.random.choice(expr.columns[1:],6,replace=False)
	for balance in tqdm.tqdm([True,False]):
		for name, geneset in zip(["target","orc-edge",'random'],[['CD274'], 
			['KIF2C','CDC25C', 'MCM6','RFC3', 'MAD2L1','SPC25'],
			random_genes]):
			
			X = np.log2(expr[geneset].values+1)
			y = resp['Response'].values

			param_grid = {'clf__C':np.arange(0.1,1,0.1)}
			
			if balance:
				cw = 'balanced'
			else:
				cw = None
			model = Pipeline([('preproc',StandardScaler()),('clf',LogisticRegression(class_weight = cw ))])

			for k in tqdm.tqdm(range(100),leave=False):
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
		plt.show()








if __name__ == '__main__':
	main({})