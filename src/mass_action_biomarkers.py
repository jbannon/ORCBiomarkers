from collections import defaultdict
import matplotlib.pyplot as plt
import argparse
from typing import Dict, Union, List
import sys
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import numpy as np 
import networkx as nx
import ot
import os
import yaml 
import utils
import pickle
import NetworkCurvature as nc
import statsmodels.stats.multitest as mt
from scipy.stats import mannwhitneyu, wilcoxon
import tqdm








def main(
	Config:Dict
	) -> None:
	
	drugs:List[str]
	tissue:str
	geneset:str
	
	
	rngSeed: int 
	alphas:List[float]
	doOneHop:bool
	minTPM:int 
	graphTops: List[str]
	weighting: str
	nIters:int 
	subsamplePct:float
	
	fdrThresh: float
	numGenes: int
	minTPM: int


	# Directory Path Stems
	dataDir:str 
	genesetDir:str
	networkDir: str 
	resultDir:str 

	normalizeWeights:bool 
	scalingConstant: Union[int,float]
	

	#Fields used throughout
	weightField:str
	renameField:str
	densityField:str 
	edgeCurvatureField:str
	nodeCurvField:str 
	normNodeCurvField:str

	sinkhorn:bool
	regParam:float
	
	
	drugs, rngSeed, alpha, doOneHop, fdrThresh, cutoff,  minTPM, geneset, nIters, subsamplePct, probThresh = \
		utils.unpack_parameters(config['EXPERIMENT_PARAMS']) 
	
	
	dataDir, genesetDir, networkDir,resultDir = utils.unpack_parameters(config['DIRECTORIES'])
	
	hopString = "OneHop" if doOneHop else "NoHop"
	
	sinkhorn, regParam = utils.unpack_parameters(config["OT_PARAMS"])
	
	normalizeWeights, scalingConstant, graphTops, weighting =\
		utils.unpack_parameters(config['NETWORK_PARAMS'])
	
	weightField, renameField, densityField, edgeCurvatureField, nodeCurvField, normNodeCurvField = \
		utils.unpack_parameters(config['FIELD_NAMES'])

	
	rng = np.random.default_rng(rngSeed)
	
	for drug in (drugProg:=tqdm.tqdm(drugs)):
		drugProg.set_description("Working on {d}".format(d=drug))
		datasetString = utils.DRUG_DATASET_MAP[drug]
		for tissue in (tissueProg := tqdm.tqdm(utils.DRUG_TISSUE_MAP[drug],leave = False)):
			tissueProg.set_description("Working on {t}".format(t=tissue))

			if drug == 'Nivo' and tissue == "KIRC":
				continue

			exprFile = utils.make_file_path(dataDir,[datasetString,drug, tissue],'expression','.csv')
			respFile = utils.make_file_path(dataDir,[datasetString,drug,tissue],'response','.csv')
			resp = pd.read_csv(respFile)
			
			expr = pd.read_csv(exprFile)
			
			trimmedExpr = expr[expr.columns[1:]]
			trimmedExpr = (trimmedExpr>=minTPM).all(axis=0)
			
			keep_genes = list(trimmedExpr[trimmedExpr].index)
			keep_cols = ['Run_ID'] + keep_genes

			expr = expr[keep_cols]

	
			gene_list, qv = utils.fetch_geneset(geneset,genesetDir,drug,tissue,fdrThresh,cutoff,probThresh)
			
			
			common_genes = [x for x in gene_list if x in expr.columns[1:]]
			
			for topology in (topProb:=tqdm.tqdm(graphTops, leave=False)):
				topProb.set_description("working on {t} topology".format(t=topology))
		
				networkFile = utils.make_file_path(networkDir,[datasetString,topology],weighting,".pickle")
		
				resPath = "".join([resultDir,"/".join(["biomarkers","mass_action",geneset,drug,tissue, str(cutoff),topology,hopString]),"/"])
				os.makedirs(resPath,exist_ok = True)
				


				with open(networkFile,"rb") as istream:
					PPI_Graph = pickle.load(istream)
	
			

				issues = list(set(expr.columns[1:]).difference(set(gene_list)))
				


				if doOneHop:
					seeds = [x for x in common_genes if x in PPI_Graph.nodes()]
					one_hop = [x for x in seeds]
					for n in seeds:
						nbrs = PPI_Graph.neighbors(n)
						one_hop.extend([x for x in nbrs if x in keep_genes])
					common_genes = [x for x in one_hop]
					# print(common_genes)

				
				LCC_Graph = utils.harmonize_graph_and_geneset(PPI_Graph,common_genes)
				
				subsetExpression = expr[['Run_ID']+[x for x in LCC_Graph.nodes()]]
			
				data = subsetExpression.merge(resp, on = 'Run_ID')


				LCC_Graph, gene_to_idx, idx_to_gene = utils.remap_LCC(LCC_Graph,renameField)
				print(LCC_Graph)
				continue
				if len(LCC_Graph.edges)<1:
					continue
				
				edgeCurvatures = defaultdict(list)
				nodeCurvatures = defaultdict(list)
				
				raw_node_curvatures = defaultdict(list)
			
				for idx, row in tqdm.tqdm(data.iterrows(),total=data.shape[0],leave=False):
				

					G = LCC_Graph.copy()
					patientResponse = row['Response']
				
					total_edge_weight = 0
				
					# print("assigning weights")
				
					for edge in G.edges():
						node1, node2 = edge[0], edge[1]
						gene1, gene2 = idx_to_gene[node1], idx_to_gene[node1]
						weight = np.round(np.log2(row[gene1]+1)*np.log2(row[gene2]+1),5)

						total_edge_weight += weight

						G[node1][node2][weightField] = weight
					
					if normalizeWeights:
						for edge in G.edges():
							node1, node2 = edge[0],edge[1]
							G[node1][node2][weightField] = scalingConstant*(G[node1][node2][weightField]/total_edge_weight)
				

					
					nc.assign_densities(G,alpha=alpha, weight=weightField,measure_name = densityField)
					
					D = nc.make_APSP_Matrix(G,weighted = True, weight=weightField)
					
					nc.compute_OR_curvature(G,D,densityField, edgeCurvatureField, sinkhorn, regParam)
					
					nc.compute_node_curvatures(G, weight = weightField,
						edge_curv = edgeCurvatureField, 
						node_curv = nodeCurvField,
						norm_node_curv = normNodeCurvField)
					
					for edge in G.edges(data=True):
						edgeCurvatures['Patient'].append(row['Run_ID'])
						edgeCurvatures['Gene1'].append(idx_to_gene[edge[0]])
						edgeCurvatures['Gene2'].append(idx_to_gene[edge[1]])
						edgeString= "{a};{b}".format(a=idx_to_gene[edge[0]] ,b = idx_to_gene[edge[1]])				
						edgeCurvatures['Edge'].append(edgeString)
						edgeCurvatures['Curvature'].append(edge[2][edgeCurvatureField])
						edgeCurvatures['Response'].append(patientResponse)

					for node in G.nodes(data= True):
						# print(node)
						nodeCurvatures['patient'].append(row['Run_ID'])
						nodeCurvatures['Gene'].append(node[1]['Gene'])
						nodeCurvatures['Raw'].append(node[1][nodeCurvField])
						nodeCurvatures['Normalized'].append(node[1][normNodeCurvField])
						nodeCurvatures['Response'].append(patientResponse)

					
				edge_data = pd.DataFrame(edgeCurvatures)
				node_data = pd.DataFrame(nodeCurvatures)
				
				edge_res = defaultdict(list)
				pvals = []
				
				for e in pd.unique(edge_data['Edge']):
					tempR = edge_data[edge_data['Edge'] == e]
					tempR = tempR[tempR['Response'] == 1]
					tempNR = edge_data[edge_data['Edge'] == e]
					tempNR = tempNR[tempNR['Response'] == 0 ]
					u, p = mannwhitneyu(tempR['Curvature'].values,tempNR['Curvature'].values,alternative = 'two-sided')
					pvals.append(p)
					edge_res['Edge'].append(e)
					edge_res['U'].append(u)
					edge_res['pval'].append(p)
			
				

				res = mt.multipletests(pvals,method = 'fdr_bh')
				edge_res['adj_pvals'] = res[1]
				edge_res = pd.DataFrame(edge_res)

				node_res = defaultdict(list)
				norm_res = defaultdict(list)
				
				node_pvals = []
				norm_pvals = []
				
				for n in pd.unique(node_data['Gene']):
					tempR = node_data[node_data['Gene'] == n]
					tempR = tempR[tempR['Response'] == 1]
					tempNR = node_data[node_data['Gene']==n]
					tempNR = tempNR[tempNR['Response']==0]
					u, p = mannwhitneyu(tempR['Normalized'].values,tempNR['Normalized'].values,alternative = 'two-sided')
					norm_pvals.append(p)
					norm_res['Gene'].append(n)
					norm_res['U'].append(u)
					norm_res['pval'].append(p)

					u, p = mannwhitneyu(tempR['Raw'].values,tempNR['Raw'].values,alternative = 'two-sided')
					node_pvals.append(p)
					node_res['Gene'].append(n)
					node_res['U'].append(u)
					node_res['pval'].append(p)

			
				res = mt.multipletests(node_pvals,method = 'fdr_bh')
				node_res['adj_pvals'] = res[1]
			
				res = mt.multipletests(norm_pvals,method = 'fdr_bh')
				norm_res['adj_pvals'] = res[1]
			
				node_res = pd.DataFrame(node_res)
				norm_res = pd.DataFrame(norm_res)
			
				statString = "alpha:\t{a}\nfdrThresh:\t{t}\nCutoff:\t{g}".format(a=alpha,t=fdrThresh,g=cutoff)
				statString = statString + "\nMinimum TPM:\t{m}".format(m=minTPM)
				statString = statString + "\nSinkhorn:\t{s}".format(s=sinkhorn)
				statString = statString + "\nEpsilon\t{e}".format(e=regParam)
				statname = resPath + "stats.txt"
				with open(statname,"w") as ostream:
					ostream.write(statString)
				
				norm_res.to_csv(resPath+"normalized_node_curvature_pvals.csv")
				
				node_res.to_csv(resPath+"node_curvature_pvals.csv")
				
				edge_res.to_csv(resPath+"edge_curvature_pvals.csv")
				
				edge_data.to_csv(resPath+"edge_curavtures.csv")
				
				node_data.to_csv(resPath+"node_curavtures.csv")
			


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)
	
	