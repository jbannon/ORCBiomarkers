import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import os 
import sys
from itertools import product


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
	fname = "../results/transfer/{d}/loo/{s}_to_{t}.csv".format(d=drug, s=s,t=t)
	
	fig_dir = "../figs/{d}/loo/{s}_to_{t}/".format(d=drug,s=s,t=t)
	df = pd.read_csv(fname)
	os.makedirs(fig_dir,exist_ok = True)
	
	for k in name2feat.keys():
		feats = name2feat[k]
		feat_df = df[df['feature'].isin(feats)]
		for pv in pd.unique(feat_df['pval']):
			p_df = feat_df[feat_df['pval']==pv]
			print(pv)
			for feat in pd.unique(p_df['feature']):
				temp = p_df[p_df['feature']==feat]
				fpr, tpr, thresholds = roc_curve(temp['true_class'].values, temp['predicted_prob'].values, pos_label=1)
				auc = roc_auc_score(temp['true_class'], temp['predicted_prob'])
				plt.plot(fpr,tpr,label = "{f} - {s}".format(f=feat,s=np.round(auc,2)))
			plt.plot([0, 1], ls="--")
			plt.legend(loc="upper left")
			plt.title("p value = {p}".format(p=pv))
			plt.savefig("{f}/{k}_{p}.png".format(f=fig_dir,k=k,p=pv))
			plt.close()
		
		