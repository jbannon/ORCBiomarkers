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



bins = 10**-2*np.arange(0,101,step = 5)
ticks = 10**-2*np.arange(0,101,step = 5)

print(["."+str(np.round(t,2)).split(".")[1] if t < 1 else str(np.round(t,2)).split(".")[0] for t in ticks ])


for drug in ["Atezo"]:
	for folder in ["full_graph","lcc_only"]:
		for tissue in ["BLCA","KIRC"]:
			for n, fn in zip(["edge","node","norm-node"], 
				["edge_curvature_pvals.csv","node_curvature_pvals.csv","normalized_node_curvature_pvals.csv"]):
				fname = "../results/biomarkers/{d}/{f}/{t}/{ff}".format(d=drug,f=folder,t=tissue,ff=fn)
				df = pd.read_csv(fname, index_col = 0)
				fig = sns.histplot(data=df, x="adj_pvals",bins = bins)
				fig.set_xticks(ticks)
				fig.set_xticklabels(["."+str(np.round(t,2)).split(".")[1] if t < 1 else str(np.round(t,2)).split(".")[0] for t in ticks ])
				plt.title("{d} {t}  {f} {n}".format(d=drug,t=tissue,f=folder,n=n))
				plt.savefig("../figs/{d}_{t}_{f}_{n}.png".format(d=drug,t=tissue, f=folder,n=n))
				plt.close()

for drug in ["Atezo"]:
	for folder in ["full_graph","lcc_only"]:
		for tissue in ["BLCA","KIRC"]:
			for n, fn in zip(["edge","node","norm-node"], 
				["edge_curvature_pvals.csv","node_curvature_pvals.csv","normalized_node_curvature_pvals.csv"]):
				fname = "../results/biomarkers/{d}/{f}/{t}/{ff}".format(d=drug,f=folder,t=tissue,ff=fn)
				df = pd.read_csv(fname, index_col = 0)
				fig = sns.histplot(data=df, x="adj_pvals",bins = bins)
				fig.set_xticks(ticks)
				fig.set_xticklabels(["."+str(np.round(t,2)).split(".")[1] if t < 1 else str(np.round(t,2)).split(".")[0] for t in ticks ])
				plt.title("{d} {t}  {f} {n}".format(d=drug,t=tissue,f=folder,n=n))
				plt.savefig("../figs/{d}_{t}_{f}_{n}.png".format(d=drug,t=tissue, f=folder,n=n))
				plt.close()


