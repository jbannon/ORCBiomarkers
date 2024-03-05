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

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2


import itertools as it
from sklearn.model_selection import LeaveOneOut
import seaborn as sns



#for drug in ['Pembro','Atezo','Nivo','Ipi+Pembro']:#,'erlotinib','crizotinib','sorafenib','sunitinib']:
for drug in ['Ipi+Pembro']:
	for graph_type in ['full_graph']:
		tissues = os.listdir("../results/biomarkers/{d}/{gt}".format(d=drug,gt=graph_type))
		for tissue in tissues:
			# if tissue!= 'KIRC':
			# 	continue
			if tissue[0]==".":
				continue
			print("\n---------------\n")
			print("\nFor Drug {d} in tissue {t}:\n".format(d=drug, t=tissue))
		
			for curve, ctype in zip(['Edge','Node','Norm Node'],['edge_curvature_pvals.csv','node_curvature_pvals.csv','normalized_node_curvature_pvals.csv']):
				if curve!= 'Norm Node':
					continue
				pvalue_file = "../results/biomarkers/{d}/{gt}/{t}/{c}".format(d=drug,gt=graph_type,t=tissue, c= ctype)
				
				pvalues = pd.read_csv(pvalue_file,index_col = 0)
				print("\n")
				print("\t{c}\n".format(c=curve))
				all_genes = []

				for thresh in [0.005, 0.01, 0.05, 0.1]:
					
					p_  = pvalues[pvalues['adj_pvals']<=thresh]
					genes = []
					

					for idx, row in p_.iterrows():
						if curve == 'Edge':
							genes.extend(row['Edge'].split(";"))
						else:
							genes.append(row['Gene'])

					genes = list(pd.unique(genes))
					new_genes = list(set(genes).difference(set(all_genes)))
					for x in new_genes:
						print(x)

					all_genes.extend(genes)
				

					print("At level {t} have {n} total genes".format(t=thresh,n=len(genes)))
					print("have added:\n")
					print(new_genes)
					print(len(new_genes))
					print("\n")
			
			sys.exit(1)	
			
	
					



