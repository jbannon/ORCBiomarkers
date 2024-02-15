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







def main():
	metrics = ["Test ROC_AUC","Test Accuracy"]
	bpath = "../results/classification/"
	for drug in os.listdir(bpath):
		# print(drug)
		dpath = bpath +drug
		if os.path.isdir(dpath):
			gpath = dpath + "/"+ "target_spread"
			if os.path.isdir(gpath):
				for tissue in os.listdir(gpath):
					tpath = gpath + "/"+tissue	
					if os.path.isdir(tpath):
						data_file = tpath + "/monte_carlo.csv"
						df = pd.read_csv(data_file,index_col = 0)
						figpath = "../figs/classification/{d}/{t}/target_spread/".format(d=drug,t=tissue)
						os.makedirs(figpath,exist_ok=True)
						for metric in metrics:
							subdf = df[df['feature'].isin(['targets','node','norm-node'])][['feature',metric]]
							order = subdf.groupby('feature').median().sort_values(by = metric).index
							fig = sns.boxplot(x="feature", y=metric,hue="feature",order = order,data=subdf, dodge = False)
							fig.set(title = "{d} {t} Performance".format(d=drug.title(),t=tissue))
							fig.set(xlabel = "Feature", ylabel = metric)
							
							metstring = metric.split(" ")[1]
							figname = figpath +metstring+".png"
							plt.savefig(figname)
							plt.close()




if __name__ == '__main__':
	main()