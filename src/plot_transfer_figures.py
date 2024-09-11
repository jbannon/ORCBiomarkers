import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import os 
import sys
from itertools import product
import scipy.stats as st

from scipy.stats import mannwhitneyu

drug = 'Atezo'
srcs = ['BLCA','KIRC']
tgts = ['BLCA','KIRC']



name2feat = {
	'edge':['edge','edge-matched-random'],
	'node':['node','node-matched-random'],
	'norm-node':['targets','norm-node','norm-node-matched-random'],
	'curvatures':['edge','node','norm-node'],
	'all':['edge','node','norm-node','KEGG','intersection']
	}


name2labels = {
	'edge':'Ricci',
	'edge-matched-random': 'Matched Random',
	'node':'Scalar',
	'norm-node':'Norm. Scalar',
	'node-matched-random':'Matched Random',
	'norm-node-matched-random': 'Matched Random',
	'KEGG':'KEGG',
	'intersection':'intersection'

	}
for s, t in product(srcs,tgts):
	if s == t:
		continue
	fname = "../results/transfer/{d}/monte_carlo/{s}_to_{t}.csv".format(d=drug, s=s,t=t)
	df = pd.read_csv(fname)
	
	
	print("{d}\n\t{s} to {t}".format(d=drug,s=s,t=t))
	# for k in ['edge']:
	for k in name2feat.keys():
		feats = name2feat[k]
		print(feats)
		feat_df = df[df['feature'].isin(feats)]
		for metric in ['Test ROC_AUC','Test Accuracy']:
			fig_dir = "../figs/{d}/monte_carlo/{s}_to_{t}/{m}/".format(d=drug,s=s,t=t,m=metric)
			os.makedirs(fig_dir,exist_ok = True)
			for pv in pd.unique(feat_df['pval']):
				temp = feat_df[feat_df['pval']==pv]
				print("for {k} -- {p}, {m}".format(k=k,p=pv,m=metric))
				exps = []
				for feat in feats:
					subdf = temp[temp['feature']==feat]
					subdf['feature'] = subdf['feature'].apply(lambda x: name2labels[x])
					# print(subdf['feature'])
					# sys.exit(1)
					exps.append(subdf[metric].values)
					print("mean = {m}".format(m=np.mean(subdf[metric])))
					interval = st.t.interval(0.95, len(subdf[metric].values)-1, loc=np.mean(subdf[metric]), scale=st.sem(subdf[metric].values))
					print("interval = {i}".format(i = [np.round(x,2) for x in interval]))
				U,p = mannwhitneyu(exps[0],exps[1])
				
				


				temp['feature'] = temp['feature'].apply(lambda x:name2labels[x])
				fig = sns.boxplot(x="feature", y=metric,hue="feature",data=temp, dodge = False)
				fig.set(title = "{d} {s} to {t} Performance\n Mann Whitney pvalue = {p}".format(d=drug.title(),s=s,t=t,p=p))
				fig.set(xlabel = "Feature", ylabel = metric)
				plt.legend([],[], frameon=False)
				plt.savefig("{f}/{k}_{p}.png".format(f=fig_dir,k=k,p=pv))
				plt.close()


sys.exit(1)


for s, t in product(srcs,tgts):
	if s == t:
		continue
	fname = "../results/transfer/{d}/loo/{s}_to_{t}.csv".format(d=drug, s=s,t=t)
	
	fig_dir = "../figs/{d}/loo/{s}_to_{t}/".format(d=drug,s=s,t=t)
	df = pd.read_csv(fname)
	os.makedirs(fig_dir,exist_ok = True)
	
	for k in name2feat.keys():
		print(k)
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


		
			
			
				
				
				
			
		
		