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
for drug in ['Nivo']:
	for graph_type in ['full_graph']:
		tissues = os.listdir("../results/biomarkers/{d}/{gt}".format(d=drug,gt=graph_type))
		for tissue in tissues:
			
			if tissue[0]==".":
				continue
			print("\n---------------\n")
			print("\nFor Drug {d} in tissue {t}:\n".format(d=drug, t=tissue))
			
			for curve, ctype in zip(['Edge'],['edge_curvature_pvals.csv']):
				
				pvalue_file = "../results/biomarkers/{d}/{gt}/{t}/{c}".format(d=drug,gt=graph_type,t=tissue, c= ctype)
				
				pvalues = pd.read_csv(pvalue_file,index_col = 0)
				
				all_genes = []
				prev_p = 0
				for thresh in [0.005, 0.01, 0.05, 0.1]:
					
					p_  = pvalues[(pvalues['adj_pvals']<=thresh) & (pvalues['adj_pvals']>prev_p)]
					prev_p = thresh
					genes = []
					print("\n")
					print(thresh)
					for e in p_['Edge']:
						print("{n1}\t{n2}	".format(n1 = e.split(";")[0],n2 = e.split(";")[1]))
			
			
			
	
					



