import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
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



def make_edge_df(pvalues:pd.DataFrame,
	thresh_vals:List[float] = [0.005, 0.01, 0.05, 0.1]
	) -> pd.DataFrame:
	

	results = defaultdict(list)
	prev_p = 0
	for thresh in thresh_vals:
		pv_subset  = pvalues[(pvalues['adj_pvals']<=thresh) & (pvalues['adj_pvals']>prev_p)]
		
		for idx, row in pv_subset.iterrows():
			v1, v2 =  row['Edge'].split(";")[0],row['Edge'].split(";")[1]
			results['Vertex 1'].append(v1)
			results['Vertex 2'].append(v2)
			results['P-val'].append(thresh)
		prev_p = thresh
	results = pd.DataFrame(results)
	return results

thresh_vals:List[float] = [0.005, 0.01, 0.05, 0.1]
for drug in ['Pembro','Atezo','Nivo','Ipi+Pembro']:

	tissues = os.listdir("../results/biomarkers/{d}/".format(d=drug))
	for tissue in tissues:
		
		if tissue[0]==".":
			continue
		if tissue in ["lcc_only","full_graph"]:
			continue
		print("\n---------------\n")
		print("\nFor Drug {d} in tissue {t}:\n".format(d=drug, t=tissue))
		os.makedirs("../results/tables/{d}/{t}".format(d=drug,t=tissue),exist_ok = True)
		
		for curve, curve_type in zip(['Ricci','Scalar','Norm-Scalar'],['edge_curvature_pvals.csv','node_curvature_pvals.csv','normalized_node_curvature_pvals.csv']):
			pvalue_file = "../results/biomarkers/{d}/{t}/{c}".format(d=drug,t=tissue, c= curve_type)
			pvalues = pd.read_csv(pvalue_file,index_col = 0)

			if curve == 'Ricci':
				edge_list = make_edge_df(pvalues)
				edge_list.to_csv("../results/tables/{d}/{t}/Edges.csv".format(d=drug,t=tissue), index = False)

			results = defaultdict(list)
			prev_p = 0
			for thresh in thresh_vals:
				pv_subset  = pvalues[(pvalues['adj_pvals']<=thresh) & (pvalues['adj_pvals']>prev_p)]
				prev_p = thresh
				for idx, row in pv_subset.iterrows():
						if curve == 'Ricci':
							results['Gene'].extend(row['Edge'].split(";"))
							results['P-val'].extend([thresh]*2)
						else:
							results['Gene'].append(row['Gene'])
							results['P-val'].extend([thresh])

			df = pd.DataFrame(results)
			df.to_csv("../results/tables/{d}/{t}/{c}_genes.csv".format(d=drug, t= tissue, c= curve))

		
	
					



