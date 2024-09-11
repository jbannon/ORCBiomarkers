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




source_drug = "Atezo"
source_tissue = "BLCA"

target_drug = "Pembro"
target_tissue = "STAD"



thresh_vals:List[float] = [0.005, 0.01, 0.05, 0.1]



clin_file = f"../data/iatlas-ici-sample_info.tsv"
clin_data = pd.read_csv(clin_file,sep="\t")


expr_file = f"../data/expression/cri/{target_drug}/{target_tissue}/expression.csv"
resp_file = f"../data/expression/cri/{target_drug}/{target_tissue}/response.csv"

expression = pd.read_csv(expr_file)
response = pd.read_csv(resp_file)
response['Response'] = response['Response'].apply(lambda x: 'Responder' if x==1 else 'Non-Responder')
ctype_map = {'edge':'edge','node':'node','norm-node':'normalized_node'}

DE_file = f"../results/biomarkers/{source_drug}/{source_tissue}/DE_genes.txt"

with open(DE_file, "r") as istream:
	DE_genes = istream.readlines()
	DE_genes = [x.rstrip() for x in DE_genes]

Filtered_DE_file = f"../results/biomarkers/{source_drug}/{source_tissue}/Filtered_DE_genes.txt"

with open(Filtered_DE_file, "r") as istream:
	Filtered_DE_genes = istream.readlines()
	Filtered_DE_genes = [x.rstrip() for x in Filtered_DE_genes]

kmf = KaplanMeierFitter()

time_col = 'PFS_d'
event_col = 'PFS_e'


for curve_type in ['edge','node','norm-node']:
	pvalue_file = f"../results/biomarkers/{source_drug}/{source_tissue}/{ctype_map[curve_type]}_curvature_pvals.csv"
	fig_dir = f"../figs/transfer/{source_drug}/{source_tissue}/{target_drug}/{target_tissue}/{curve_type}/"
	os.makedirs(fig_dir,exist_ok=True)
	pvalues = pd.read_csv(pvalue_file,index_col = 0)

	for thresh in thresh_vals:
		
		pv_subset  = pvalues[(pvalues['adj_pvals']<=thresh)]
		
		if curve_type == 'edge':
			genes = [ x.split(";") for x in pv_subset['Edge']]
			genes = [gene for pairs in genes for gene in pairs]
		else:
			genes = pv_subset['Gene'].values
		
		genes = [x for x in genes if x not in Filtered_DE_genes]
		
		genes = sorted(list(set(genes)))
		genes = [x for x in genes if x in expression.columns]

		if len(genes)==0:
			continue

		expr_ = expression[['Run_ID'] + genes]
		pre = StandardScaler()
		X = np.log2(expr_[expr_.columns[1:]].values+1)
		X  = pre.fit_transform(X)

		best_score = -np.inf
		best_k = 0
	
		for k in [2,3,4]:
			
			# clustering = KMeans(n_clusters=k).fit(X)
			
			kmeans = KMeans(n_clusters=k, random_state=42)# print(thresh)
			score = silhouette_score(X, kmeans.fit_predict(X))
			if score >best_score:
				best_k = k
				best_score = score
		
		clustering = KMeans(n_clusters = best_k).fit(X)
		

		label_df = pd.DataFrame({'Run_ID':expression['Run_ID'].values,'Cluster':clustering.labels_})
		df = label_df.merge(clin_data,on='Run_ID')
		df.dropna(inplace=True,subset = [time_col,event_col])
		
		ax = plt.subplot(111)
		for cluster in sorted(pd.unique(df['Cluster'])):
			# print(cluster)
			
			kmf.fit(df[df['Cluster']==cluster][time_col], 
				event_observed=df[df['Cluster']==cluster][event_col], 
				label=f"Cluster {cluster}")
			kmf.plot_survival_function(ax=ax,ci_show=False)
			# at_risk_counts=True)

		results = logrank_test(df[df['Cluster']==0][time_col], 
			df[df['Cluster']==1][time_col], 
			df[df['Cluster']==0][event_col],
			df[df['Cluster']==1][event_col])
		
		ax.set_xlabel('Time (Days)', fontsize = 18)
		ax.set_ylabel("Survival Probability",fontsize = 18)
		ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
		ax.set_yticklabels(ax.get_yticklabels(), fontsize=15)
		ax.set_title(f"({source_drug}, {source_tissue}) to ({target_drug}, {target_tissue})", {"fontsize":20})

		legend = ax.legend()
		legend.remove()
		plt.tight_layout()
		res_pv = results.p_value
		with open(f"{fig_dir}{curve_type}_{thresh}_pvalue.txt","w") as ostream:
			ostream.writelines([str(res_pv)+"\n"])
		# kmf.plot_survival_function(ax=ax)
		# print(f"{drug} {tissue} {thresh}")
		# print(res_pv)
		# plt.show()
		# sys.exit(1)
		# plt.title(f"{source_drug} --> {target_drug}")
		plt.tight_layout()
		plt.savefig(f"{fig_dir}{curve_type}_{thresh}.png",dpi=300)
		plt.close()



