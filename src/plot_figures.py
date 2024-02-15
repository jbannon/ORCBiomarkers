import numpy as np
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
	'edge':['targets','edge','edge-random','edge+targets'],
	'node':['targets','node','node-random','node+targets'],
	'norm-node':['targets','node-norm','node-norm-random','node-norm+targets'],
	'curvatures':['targets','edge','node','node-norm'],
	'target_curvatures':['targets','edge+targets','node+targets','node-norm+targets']}
	metrics = ["Test ROC_AUC","Test Accuracy"]
	bpath = "../results/classification/"

	for drug in os.listdir(bpath):
		# print(drug)
		dpath = bpath +drug
		if os.path.isdir(dpath):
			for graph_type in os.listdir(bpath+drug):
				gpath = dpath + "/"+ graph_type
				if graph_type == "target_spread":
					continue
				if os.path.isdir(gpath):

					for tissue in os.listdir(gpath):
						tpath = gpath + "/"+tissue
						
						if os.path.isdir(tpath):
							data = tpath + "/monte_carlo.csv"
							df = pd.read_csv(data, index_col = 0)
							pre_drop_shape = df.shape[0]
							df = df[df['Test TP']>=0]
							post_drop_shape = df.shape[0]

							dim_dict = {}
							for feature_type in pd.unique(df['feature']):
								feat_df = df[df['feature']==feature_type]
								# print(feat_df.columns)
								# sys.exit(1)
								# dim_dict[feature_type] = np.mean(feat_df['feature_dim'].values)
								
								
							
							figpath = "../figs/classification/{d}/{t}/{e}/".format(d=drug,t=tissue,e=graph_type)
							drop_stats = "Before: {n}\nAfter: {nn}".format(n=pre_drop_shape,nn=post_drop_shape)
							with open(figpath+"dropstats.txt","w") as ostream:
								ostream.write(drop_stats)
							
							# with open(figpath + "dim_stats.txt","w") as ostream:
							# 	ostream.writelines([k+"\t"+str(dim_dict[k])+"\n" for k in dim_dict.keys()])
							
							for metric in metrics:
								for k in Ls.keys():
									subdf = df[df['feature'].isin(Ls[k])][['feature',metric]]
									order = subdf.groupby('feature').median().sort_values(by = metric).index
									
									fig = sns.boxplot(x="feature", y=metric,hue="feature",order = order,data=subdf, dodge = False)
									gstring = "LCC" if graph_type =='lcc_only' else "Full"

									fig.set(title = "{d} {t} {l} Performance".format(d=drug.title(),t=tissue, l = gstring))
									fig.set(xlabel = "Feature", ylabel = metric)

									
									os.makedirs(figpath,exist_ok=True)
									metstring = metric.split(" ")[1]
									
									figname = figpath + k+"_"+metstring+".png"
									plt.savefig(figname)
									plt.close()
									
									
								
										
										
										

if __name__ == '__main__':
	main()