import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os

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


import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2

import itertools as it
import seaborn as sns


#for drug in ['Pembro','Atezo','Nivo','Ipi+Pembro']:
for drug in ['Atezo']:
	tissues = os.listdir("../results/biomarkers/{d}/".format(d=drug))
	for tissue in tissues:
			
		if tissue[0]==".":
			continue
		if tissue in ["lcc_only","full_graph"]:
			continue
		print("\n---------------\n")
		print("\nFor Drug {d} in tissue {t}:\n".format(d=drug, t=tissue))
		
		results = defaultdict(list)
		for curve, curve_type in zip(['Edge'],['edge_curvature_pvals.csv']):
			pvalue_file = "../results/biomarkers/{d}/{t}/{c}".format(d=drug,t=tissue, c= curve_type)
			pvalues = pd.read_csv(pvalue_file,index_col = 0)
			all_genes = []
			prev_p = 0

			for thresh in [0.005, 0.01, 0.05, 0.1]:
				pv_subset  = pvalues[(pvalues['adj_pvals']<=thresh) & (pvalues['adj_pvals']>prev_p)]
				prev_p = thresh
				genes = []
				for idx, row in pv_subset.iterrows():
					v1, v2 =  row['Edge'].split(";")[0],row['Edge'].split(";")[1]
					results['Vertex 1'].append(v1)
					results['Vertex 2'].append(v2)
					results['P-val'].append(thresh)
					# print("{n1}\t{n2}	".format(n1 = e.split(";")[0],n2 = e.split(";")[1]))
		df = pd.DataFrame(results)
		df.to_csv("../figs/{d}_{t}_edges.csv".format(d=drug,t=tissue), index = False)
		
			
			
	
					



