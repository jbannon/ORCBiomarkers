import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import os 
import sys
from itertools import product
import scipy.stats as st

drug = 'Atezo'
srcs = ['BLCA','KIRC']
tgts = ['BLCA','KIRC']



name2feat = {
	'edge':['edge','edge-matched-random'],
	'node':['node','node-matched-random'],
	'norm-node':['targets','norm-node','norm-node-matched-random'],
	'curvatures':['edge','node','norm-node']
	}



for s, t in product(srcs,tgts):
	if s == t:
		continue
	fname = "../results/transfer/{d}/monte_carlo/{s}_to_{t}.csv".format(d=drug, s=s,t=t)
	df = pd.read_csv(fname)
	
	df = df[df['type']=='full_graph']
	
	
	print("{d}\n\t{s} to {t}".format(d=drug,s=s,t=t))
	# for k in ['edge']:
	for k in name2feat.keys():
		feats = name2feat[k]
		feat_df = df[df['feature'].isin(feats)]
		for metric in ['Test ROC_AUC','Test Accuracy']:
			fig_dir = "../figs/{d}/monte_carlo/{s}_to_{t}/{m}/".format(d=drug,s=s,t=t,m=metric)
			os.makedirs(fig_dir,exist_ok = True)
			for pv in pd.unique(feat_df['pval']):
				temp = feat_df[feat_df['pval']==pv]
				print("for {k} -- {p}, {m}".format(k=k,p=pv,m=metric))
				for feat in feats:
					subdf = temp[temp['feature']==feat]
					print(feat)
					print("mean = {m}".format(m=np.mean(subdf[metric])))
					interval = st.t.interval(0.95, len(subdf[metric].values)-1, loc=np.mean(subdf[metric]), scale=st.sem(subdf[metric].values))
					print("interval = {i}".format(i = [np.round(x,2) for x in interval]))
				print("\n")

				
				fig = sns.boxplot(x="feature", y=metric,hue="feature",data=temp, dodge = False)
				fig.set(title = "{d} {s} to {t} Performance\n pvalue = {p}".format(d=drug.title(),s=s,t=t,p=pv))
				fig.set(xlabel = "Feature", ylabel = metric)
				plt.legend([],[], frameon=False)
				plt.savefig("{f}/{k}_{p}.png".format(f=fig_dir,k=k,p=pv))
				plt.close()


		
			
			
				
				
				
			
		
		