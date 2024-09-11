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



src_drug = "Pembro"
tgt_drug = "Pembro"

src_tissue = "STAD"
tgt_tissue = "SKCM"


pvs = [0.01, 0.05, 0.1]
# pvs = [0.1]

for pv in pvs:
	res = f"../results/transfer/{src_drug}/{src_tissue}/{tgt_drug}/{tgt_tissue}/monte_carlo_{pv}.csv"
	res = pd.read_csv(res)
	ax = sns.boxplot(data=res, x='geneset',y='Test Accuracy')
	plt.title(f"{pv}")
	plt.show()
