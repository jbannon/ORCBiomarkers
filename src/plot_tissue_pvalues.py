import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
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
import seaborn as sns



#for drug in ['Pembro','Atezo','Nivo','Ipi+Pembro']:#,'erlotinib','crizotinib','sorafenib','sunitinib']:
for drug in ['Atezo']:
	
	for graph_type in ['full_graph']:
		tissues = os.listdir("../results/biomarkers/{d}/{gt}".format(d=drug,gt=graph_type))
		for tissue in tissues:
			if tissue[0]==".":
				continue
			print("\nFor Drug {d} in tissue {t}:\n".format(d=drug, t=tissue))
			for curve, ctype in zip(['Edge','Node','Norm Node'],['edge_curvature_pvals.csv','node_curvature_pvals.csv','normalized_node_curvature_pvals.csv']):
				if curve!='Norm Node':
					continue
				pvalue_file = "../results/biomarkers/{d}/{gt}/{t}/{c}".format(d=drug,gt=graph_type,t=tissue, c= ctype)
				pvalues = pd.read_csv(pvalue_file,index_col = 0)
				print("\n")
				print("\t{c}\n".format(c=curve))
				
				for thresh in [0.01, 0.005]:
					p_  = pvalues[pvalues['adj_pvals']<=thresh]
					genes = []
					print(p_.sort_values('adj_pvals'))
					for idx, row in p_.iterrows():
						if curve == 'Edge':
							genes.extend(row['Edge'].split(";"))
						else:
							genes.append(row['Gene'])

					genes = list(pd.unique(genes))
					print(genes)
					print("{t}\t{n}".format(t=thresh,n=len(genes)))
				
			
	print(["*"]*20)
					



