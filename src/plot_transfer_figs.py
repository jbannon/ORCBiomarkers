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
	'edge':['edge','edge-matched-random'],
	'node':['node','node-matched-random'],
	'norm-node':['targets','norm-node','norm-node-matched-random'],
	'curvatures':['edge','node','norm-node']}

	metrics = ["Test ROC_AUC","Test Accuracy"]
	exps = {'Erlotinib (BRCA->LUAD) Transfer':("../results/transfer/erlotinib/BRCA_to_LUAD.csv",'Erlotinib'),
			'Atezo (BLCA->KIRC) Transfer': ("../results/transfer/Atezo/BLCA_to_KIRC.csv",'Atezo')
	}

	for key in exps.keys():
		df = pd.read_csv(exps[key][0],index_col = 0)
		for et in ['strict','loose']:
			for metric in metrics:

				for pval in [0.1] if key =="Erlotinib (BRCA->LUAD) Transfer'" else [0.1,0.05]:
					for k in Ls.keys():
						subdf = df[(df['pval']==pval) &(df['type']==et)]
						subdf = subdf[subdf['feature'].isin(Ls[k])]
						if subdf.shape[0]==0:
							continue
						order = subdf.groupby('feature').median().sort_values(by = metric).index
						fig = sns.boxplot(x="feature", y=metric, hue="feature",data=subdf,order = order, dodge = False)
						fig.set(title = key + " "+str(pval))
						plt.legend([],[], frameon=False)
						fig.set(xlabel = "Feature")
						fig.set(ylabel = metric)
						plt.savefig("../figs/{k}_{e}_{p}_{kk}_{m}_transfer.png".format(k=key, e = et, p = pval, kk = k,m=metric))
						plt.close()
						
										
										
										

if __name__ == '__main__':
	main()