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



from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore")
import os
from collections import defaultdict
from typing import Dict, Union, List
import sys
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import numpy as np 
from matplotlib_venn import venn3, venn2
import itertools as it
import seaborn as sns
from sklearn.cluster import SpectralClustering,KMeans
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test




# drug = 'Atezo'
# tissues = ['BLCA','KIRC']
# thresh = 0.005


# drug = 'Pembro'
# tissues = ['STAD','SKCM']
# thresh = 0.01


drug = 'Nivo'
tissues = ['SKCM']
thresh = 0.05

# curve_type = 'edge'
# union_genes = ['CD274','PDCD1','CTLA4'] # start with drug targets
union_genes = []
for tissue in tissues:
	
	for curve_type in ['edge','node','normalized_node']:
		
		
		DE_file = f"../results/biomarkers/{drug}/{tissue}/DE_genes.txt"

		with open(DE_file, "r") as istream:
			DE_genes = istream.readlines()
		DE_genes = [x.rstrip() for x in DE_genes]
	
	

		pvalue_file = f"../results/biomarkers/{drug}/{tissue}/{curve_type}_curvature_pvals.csv"

		pvalues = pd.read_csv(pvalue_file,index_col = 0)
		
		pv_subset  = pvalues[(pvalues['adj_pvals']<=thresh)]

		if curve_type == 'edge':
			genes = [ x.split(";") for x in pv_subset['Edge']]
			genes = [gene for pairs in genes for gene in pairs]
		else:
			genes = pv_subset['Gene'].values

		


		genes = [x for x in genes if x not in DE_genes]
		union_genes.extend(genes)

union_genes = list(pd.unique(union_genes))

votes = {}
drug = 'Nivo'
tissues = ['SKCM']
thresh = 0.05
for tissue in tissues:
	expression_file = f"../data/expression/cri/{drug}/{tissue}/expression.csv"
	response_file = f"../data/expression/cri/{drug}/{tissue}/response.csv"
	expression = pd.read_csv(expression_file)
	response = pd.read_csv(response_file)
	

	for g in union_genes:
		expression[g]=np.log2(expression[g]+1)

	y = response['Response'].values
	X = expression[union_genes]

	for name, criterion in zip(['mut-info','chi2','f'],[mutual_info_classif,chi2,f_classif]):
		selector = SelectKBest(mutual_info_classif, k = 'all') 
		key = f"{tissue} - {name}"
		X_new = selector.fit_transform(X,y)
		indices = np.argsort(-selector.scores_)

		ranked_list = [union_genes[i] for i in indices]
		votes[key] = ranked_list



borda_counts = {x:[] for x in union_genes}
for voter in votes.keys():
	rank_list = votes[voter]
	for j in range(len(rank_list)):
		gene = ranked_list[j]
		borda_counts[gene].append(len(rank_list)-j)

summed_counts = {x:np.sum(borda_counts[x]) for x in borda_counts.keys()}
sorted_counts = {k: v for k, v in sorted(summed_counts.items(), key=lambda item: item[1])}
aggregate_list = list(sorted_counts.keys())
aggregate_list.reverse()

with open(f"../results/biomarkers/{drug}/aggregate_list.txt","w") as ostream:
	ostream.writelines([x+"\n" for x in aggregate_list])

		


