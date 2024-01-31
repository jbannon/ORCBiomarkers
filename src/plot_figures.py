import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from typing import Dict, Union, List
import pandas as pd
import os
import yaml 
import utils
import pickle









def main(
	# Config:Dict
	) -> None:
	

	Ls = {
	'edge':['targets','edge','edge-random'],
	'node':['targets','node','node-random'],
	'norm-node':['targets','node-norm','node-norm-random'],
	'curvatures':['targets','edge','node','node-norm']}
	metrics = ["Test ROC_AUC","Test Accuracy"]
	bpath = "../results/classification/"

	for drug in os.listdir(bpath):
		# print(drug)
		dpath = bpath +drug
		if os.path.isdir(dpath):
			for tissue in os.listdir(bpath+drug):
				tpath = dpath + "/"+tissue
				if os.path.isdir(tpath):
					for etype in os.listdir(tpath):
						epath = tpath +"/" + etype
						if os.path.isdir(epath):
							for f in os.listdir(epath):
								fname = epath+"/"+f
								df = pd.read_csv(fname,index_col = 0)
								df = df[df['feature_dim']>=0]
								for metric in metrics:
									for k in Ls.keys():
										subdf = df[df['feature'].isin(Ls[k])][['feature',metric]]
										order = subdf.groupby('feature').median().sort_values(by = metric).index
										
										fig = sns.boxplot(x="feature", y=metric,hue="feature",order = order,data=subdf, dodge = False)

										fig.set(title = "{d} {t} {l} Performance".format(d=drug.title(),t=tissue, l = etype.title()))
										fig.set(xlabel = "Feature", ylabel = metric)
										figpath = "../figs/classification/{d}/{t}/{e}/".format(d=drug,t=tissue,e=etype)
										os.makedirs(figpath,exist_ok=True)
										figname = figpath + k+"_"+metric+".png"
										plt.savefig(figname)
										plt.close()
										
										
										

if __name__ == '__main__':
	main()