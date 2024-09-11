import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
# warnings.simplefilter(action='ignore', category=SettingWithCopyWarning)
import os
from collections import defaultdict
from typing import Dict, Union, List
import sys
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import numpy as np 
from matplotlib_venn import venn3, venn2
import itertools as it
import seaborn as sns
import statsmodels.stats.multitest as mt
from scipy.stats import mannwhitneyu, wilcoxon


drug = 'Nivo'
tissue = 'SKCM'



thresh_vals:List[float] = [0.005, 0.01, 0.05, 0.1]


expr_file = f"../data/expression/cri/{drug}/{tissue}/expression.csv"
resp_file = f"../data/expression/cri/{drug}/{tissue}/response.csv"

expression = pd.read_csv(expr_file)
response = pd.read_csv(resp_file)
response['Response'] = response['Response'].apply(lambda x: 'Responder' if x==1 else 'Non-Responder')

DE_file = f"../results/biomarkers/{drug}/{tissue}/DE_genes.txt"

with open(DE_file, "r") as istream:
	DE_genes = istream.readlines()
	DE_genes = [x.rstrip() for x in DE_genes]

Filtered_DE_file = f"../results/biomarkers/{drug}/{tissue}/Filtered_DE_genes.txt"

with open(Filtered_DE_file, "r") as istream:
	Filtered_DE_genes = istream.readlines()
	Filtered_DE_genes = [x.rstrip() for x in Filtered_DE_genes]
print(DE_genes)
print(Filtered_DE_genes)


for curve_type in ["edge","node","normalized_node"]:
	

	pvalue_file = f"../results/biomarkers/{drug}/{tissue}/{curve_type}_curvature_pvals.csv"
	fig_dir = f"../figs/{drug}/{tissue}/{curve_type}/"
	res_dir = f"../results/biomarkers/{drug}/{tissue}/{curve_type}/"
	os.makedirs(fig_dir, exist_ok=True)
	os.makedirs(res_dir, exist_ok=True)

	pvalues = pd.read_csv(pvalue_file,index_col = 0)

	print("*********")
	print(fig_dir)
	print("*********")

	
	
	prev_list = []
	for thresh in thresh_vals:
		

		pv_subset  = pvalues[(pvalues['adj_pvals']<=thresh)]
		
		if curve_type == 'edge':
			genes = [ x.split(";") for x in pv_subset['Edge']]
			genes = [gene for pairs in genes for gene in pairs]
		else:
			genes = pv_subset['Gene'].values
		
		genes = sorted(list(set(genes)))
		genes = [x for x in genes if x in expression.columns]
		
		expr_ = expression[['Run_ID']+genes]	
		for g in genes:
			expr_[g]=np.log2(expr_[g]+1)

		

		data = expr_.merge(response, on='Run_ID')
		
		if len(genes)==0:
			continue
			
		data.drop(columns = ['Run_ID'],inplace=True)
		data = data.melt(id_vars=['Response'],var_name = 'Gene')

		genes_in_seed_genes = [x for x in list(set(data['Gene'].values)) if x in DE_genes]
		genes_in_kept_seed_genes = [x for x in list(set(data['Gene'].values)) if x in Filtered_DE_genes]


		with open(f"{res_dir}genes_in_seed_genes_{thresh}.txt","w") as ostream:
			ostream.writelines([str(x)+"\n" for x in genes_in_seed_genes])

		with open(f"{res_dir}genes_in_kept_seed_genes_{thresh}.txt","w") as ostream:
			ostream.writelines([str(x)+"\n" for x in genes_in_seed_genes])

		with open(f"{res_dir}genes_{thresh}.txt","w") as ostream:
			ostream.writelines([str(x)+"\n" for x in genes])

		p_vals = []
		adj_pvals = []
		_g = []
		
		for gene in pd.unique(data['Gene']):
			temp = data[data['Gene']==gene]
			temp_r =temp[temp['Response']=='Responder']
			temp_nr =temp[temp['Response']=='Non-Responder']
			r_exp = temp_r['value'].values
			nr_exp = temp_nr['value'].values
			U,p = mannwhitneyu(r_exp, nr_exp)
			p_vals.append(p)
			_g.append(gene)

		adj_pvals = mt.multipletests(p_vals,method = 'fdr_bh')
		
		res_string = []
		for _t in [0.005,0.01,0.05,0.1]:
			_kept_genes = [_g[i] for i in np.where(adj_pvals[1]<=_t)[0]]
			
			res_string.append(f"At threshold {_t} we have the genes: {[str(x)+' ' for x in _kept_genes]}")
			res_string.append("\n")
			_kept_genes = [x for x in _kept_genes if x not in DE_genes]
			res_string.append(f"At threshold {_t}  not in DE genes: {[str(x)+' ' for x in _kept_genes]}")
			res_string.append("\n")
			_kept_genes = [x for x in _kept_genes if x not in Filtered_DE_genes]
			res_string.append(f"At threshold {_t}  not in Filtered DE genes: {[str(x)+' ' for x in _kept_genes]}")
			res_string.append("\n")
		

		new_genes = [x for x in _kept_genes if x not in prev_list]

		with open(f"{res_dir}new_genes_at_{thresh}.txt","w") as ostream:
			ostream.writelines([str(x)+"\n" for x in new_genes])
		prev_list = [x for x in _kept_genes]
		

		with open(f"{res_dir}biomarkers_{thresh}.txt","w") as ostream:
			ostream.writelines(res_string)
			

		# ---
		# ax.set_title(f"{drug.title()}  {tissue}  {model_name_map[model]}", fontsize = 18)
		
		
		
				# sns.set(font_scale = 1.3)
				# ax.set(title = f"{drug.title()}  {tissue}  {model_name_map[model]}",
				#  	xlabel = 'Covariates',
				# 	ylabel = metric_map[metric])
				
				# plt.tight_layout()
				#palette = {'Responder':'orange','Non-Responder':'blue'})
		ax = sns.boxplot(data=data,y='value', x='Gene',hue= 'Response')
		ax.set_xlabel('Gene', fontsize = 18)
		ax.set_ylabel("Log2 TPM+1 Expression",fontsize = 18)
		ax.set_xticklabels(ax.get_xticklabels(), rotation=35,fontsize=12,ha="right")
		# ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
		ax.set_yticklabels(ax.get_yticklabels(), fontsize=15)
		# ax.tick_params('x',labelsize=15,labelrotation=15)
		# ax.tick_params('y',labelsize=15)
		legend = ax.legend()
		legend.remove()
		# plt.legend(ncol=2,frameon=False)

		plt.tight_layout()
		fname = f"{fig_dir}curvature_{thresh}.png"
		plt.savefig(fname,dpi=300)
		# plt.show()
		plt.close()
		

		data_sub = data[data['Gene'].isin(_kept_genes)]
		ax = sns.boxplot(data=data_sub,y='value', x='Gene',hue= 'Response')
		ax.set_xlabel('Gene', fontsize = 18)
		ax.set_ylabel("Log2 TPM+1 Expression",fontsize = 18)
		ax.set_xticklabels(ax.get_xticklabels(), rotation=35,fontsize=12,ha="right")
		# ax.tick_params('x',labelsize=15,labelrotation=45)
		ax.tick_params('y',labelsize=12)
		legend = ax.legend()
		legend.remove()

		plt.tight_layout()
		fname = f"{fig_dir}curvature_{thresh}_sub.png"
		plt.savefig(fname,dpi=300)
		plt.close()


# sys.exit(1)
# curve_type = 'node_curvature_pvals.csv'
# pvalue_file = "../results/biomarkers/{d}/{t}/{c}".format(d=drug,t=tissue, c= curve_type)
# thresh_vals:List[float] = [0.005, 0.01, 0.05, 0.1]	


# expr_file = f"../data/expression/cri/{drug}/{tissue}/expression.csv"
# resp_file = f"../data/expression/cri/{drug}/{tissue}/response.csv"

# expression = pd.read_csv(expr_file)
# response = pd.read_csv(resp_file)
# pvalues = pd.read_csv(pvalue_file,index_col = 0)
# response['Response'] = response['Response'].apply(lambda x: 'Responder' if x==1 else 'Non-Responder')

# for thresh in thresh_vals:
# 	# print(thresh)
# 	pv_subset  = pvalues[(pvalues['adj_pvals']<=thresh)]
# 	# print(pv_subset)
# 	# genes = [ x.split(";") for x in pv_subset['Edge']]
# 	genes = pv_subset['Gene'].values
# 	# genes = [gene for pairs in genes for gene in pairs]
# 	genes = sorted(list(set(genes)))
# 	genes = [x for x in genes if x in expression.columns]
# 	if len(genes)==0:
# 		continue
# 	expr_ = expression[['Run_ID']+genes]	
# 	for g in genes:
# 		expr_[g]=np.log2(expr_[g]+1)

	

# 	data = expr_.merge(response, on='Run_ID')

# 	data.drop(columns = ['Run_ID'],inplace=True)

# 	data = data.melt(id_vars=['Response'],var_name = 'Gene')
# 	# print(data['Response'].value_counts())
# 	# data = data[data['Gene'].isin(['CXCL9','CXCL11','CCRI','CD8A','GZMB'])]
# 	p_vals = []
# 	adj_pvals = []
# 	_g = []
# 	for gene in pd.unique(data['Gene']):
# 		temp = data[data['Gene']==gene]
# 		temp_r =temp[temp['Response']=='Responder']
# 		temp_nr =temp[temp['Response']=='Non-Responder']
# 		r_exp = temp_r['value'].values
# 		nr_exp = temp_nr['value'].values
# 		U,p = mannwhitneyu(r_exp, nr_exp)
# 		p_vals.append(p)
# 		_g.append(gene)

# 	adj_pvals = mt.multipletests(p_vals,method = 'fdr_bh')
	
# 	res_string = []
# 	for _t in [0.005,0.01,0.05,0.1]:
# 		_kept_genes = [_g[i] for i in np.where(adj_pvals[1]<=_t)[0]]
		
# 		res_string.append(f"At threshold {_t} we have the genes: {[str(x)+' ' for x in _kept_genes]}")
# 		res_string.append("\n")
	
# 	with open(f"../results/biomarkers/{drug}/{tissue}/node_biomarkers_{thresh}.txt","w") as ostream:
# 		ostream.writelines(res_string)


# 	ax = sns.boxplot(data=data,y='value', x='Gene',hue= 'Response')
# 	ax.set(ylabel = "Log2 TPM+1 Expression", xlabel = "Gene")
# 	plt.xticks(rotation=45)

# 	plt.legend(ncol=2,frameon=False)

# 	# sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

# 	plt.tight_layout()
# 	title = f"{fig_dir}_node_{thresh}.png"
# 	# plt.show()
# 	plt.savefig(title)
# 	plt.close()


# curve_type = 'normalized_node_curvature_pvals.csv'
# pvalue_file = "../results/biomarkers/{d}/{t}/{c}".format(d=drug,t=tissue, c= curve_type)
# thresh_vals:List[float] = [0.005, 0.01, 0.05, 0.1]


# expr_file = f"../data/expression/cri/{drug}/{tissue}/expression.csv"
# resp_file = f"../data/expression/cri/{drug}/{tissue}/response.csv"

# expression = pd.read_csv(expr_file)
# response = pd.read_csv(resp_file)
# pvalues = pd.read_csv(pvalue_file,index_col = 0)
# response['Response'] = response['Response'].apply(lambda x: 'Responder' if x==1 else 'Non-Responder')

# for thresh in thresh_vals:
# 	# print(thresh)
# 	pv_subset  = pvalues[(pvalues['adj_pvals']<=thresh)]
# 	# print(pv_subset)
# 	# genes = [ x.split(";") for x in pv_subset['Edge']]
# 	genes = pv_subset['Gene'].values
# 	# genes = [gene for pairs in genes for gene in pairs]
# 	genes = sorted(list(set(genes)))
# 	genes = [x for x in genes if x in expression.columns]
# 	if len(genes)==0:
# 		continue
# 	expr_ = expression[['Run_ID']+genes]	
# 	for g in genes:
# 		expr_[g]=np.log2(expr_[g]+1)

	

# 	data = expr_.merge(response, on='Run_ID')

# 	data.drop(columns = ['Run_ID'],inplace=True)

# 	data = data.melt(id_vars=['Response'],var_name = 'Gene')
# 	# print(data['Response'].value_counts())
# 	# data = data[data['Gene'].isin(['CXCL9','CXCL11','CCRI','CD8A','GZMB'])]
# 	p_vals = []
# 	adj_pvals = []
# 	_g = []
# 	for gene in pd.unique(data['Gene']):
# 		temp = data[data['Gene']==gene]
# 		temp_r =temp[temp['Response']=='Responder']
# 		temp_nr =temp[temp['Response']=='Non-Responder']
# 		r_exp = temp_r['value'].values
# 		nr_exp = temp_nr['value'].values
# 		U,p = mannwhitneyu(r_exp, nr_exp)
# 		p_vals.append(p)
# 		_g.append(gene)

# 	adj_pvals = mt.multipletests(p_vals,method = 'fdr_bh')
	
# 	res_string = []
# 	for _t in [0.005,0.01,0.05,0.1]:
# 		_kept_genes = [_g[i] for i in np.where(adj_pvals[1]<=_t)[0]]
		
# 		res_string.append(f"At threshold {_t} we have the genes: {[str(x)+' ' for x in _kept_genes]}")
# 		res_string.append("\n")
	
# 	with open(f"../results/biomarkers/{drug}/{tissue}/norm_node_biomarkers_{thresh}.txt","w") as ostream:
# 		ostream.writelines(res_string)

# 	ax = sns.boxplot(data=data,y='value', x='Gene',hue= 'Response')
# 	ax.set(ylabel = "Log2 TPM+1 Expression", xlabel = "Gene")
# 	plt.xticks(rotation=45)

# 	plt.legend(ncol=2,frameon=False)

# 	# sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

# 	plt.tight_layout()
# 	title = f"../figs/{drug}/{drug}_{tissue}_norm_node_{thresh}.png"
# 	# plt.show()
# 	plt.savefig(title)
# 	plt.close()

	
					



