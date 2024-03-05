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
	
	
	

	metrics = ["Test ROC_AUC","Test Accuracy"]
	exps = {
			'Atezo (BLCA->KIRC) Transfer': ("../results/transfer/Atezo/BLCA_to_KIRC.csv",'Atezo')
	}
	for key in exps.keys():
		df = pd.read_csv(exps[key][0],index_col = 0)
		for et in ['full_graph','lcc_only']:
			for metric in metrics:
				for pval in[0.1,0.05,0.01,0.005]:
					for k in Ls.keys():
						subdf = df[(df['pval']==pval) &(df['type']==et)]
						subdf = subdf[subdf['feature'].isin(Ls[k])]
						if subdf.shape[0]==0:
							continue
						fig = sns.boxplot(x="feature", y=metric, hue="feature",data=subdf, dodge = False)
						fig.set(title = key + " "+str(pval))
						plt.legend([],[], frameon=False)
						fig.set(xlabel = "Feature")
						fig.set(ylabel = metric)
						plt.show()
						# sys.exit(1)
						# plt.savefig("../figs/{k}_{e}_{p}_{kk}_{m}_transfer.png".format(k=key, e = et, p = pval, kk = k,m=metric))
						# plt.close()
						
										
										
										

if __name__ == '__main__':
	main()