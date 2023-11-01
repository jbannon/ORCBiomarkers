from scipy import stats
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split, GridSearchCV
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
# import statsmodels.stats.multitest.multipletests
import pickle
import NetworkCurvature as nc
import statsmodels.stats.multitest as mt
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.linear_model import LogisticRegression
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
	
	
	
	drugs, rngSeed, alpha, doOneHop, fdrThresh, numGenes,  minTPM, geneset, nIters, subsamplePct = \
		utils.unpack_parameters(config['EXPERIMENT_PARAMS']) 
	
	dataDir, genesetDir, networkDir,resultDir = utils.unpack_parameters(config['DIRECTORIES'])
	

	sinkhorn, regParam = utils.unpack_parameters(config["OT_PARAMS"])
	normalizeWeights, scalingConstant, graphTops, weighting =\
		utils.unpack_parameters(config['NETWORK_PARAMS'])
	
	weightField, renameField, densityField, edgeCurvatureField, nodeCurvField, normNodeCurvField = \
		utils.unpack_parameters(config['FIELD_NAMES'])

	
	
	# np.random.seed(rngSeed)

	
	
	

	
	rng = np.random.default_rng(rngSeed)
	# for drug in (drugProg:=tqdm.tqdm(drugs)):
	for drug in drugs:
		print(drug)
		# drugProg.set_description("Working on {d}".format(d=drug))
		
		# for tissue in (tissueProg := tqdm.tqdm(utils.DRUG_TISSUE_MAP[drug],leave = False)):
		for tissue in utils.DRUG_TISSUE_MAP[drug]:
			# tissueProg.set_description("Working on {t}".format(t=tissue))
			if drug == 'Nivo' and tissue == 'KIRC' and geneset =='DE':
				continue

			exprFile = utils.make_file_path(dataDir,[drug, tissue],'expression','.csv')
			# print(exprFile)
			respFile = utils.make_file_path(dataDir,[drug,tissue],'response','.csv')
			
			# print(diffExpFile)
			# print(respFile)


			resp = pd.read_csv(respFile)
			
			
			expr = pd.read_csv(exprFile)
			# print(expr.shape)
			trimmedExpr = expr[expr.columns[1:]]
			trimmedExpr = (trimmedExpr>=minTPM).all(axis=0)
			keep_genes = list(trimmedExpr[trimmedExpr].index)
			keep_cols = ['Run_ID'] + keep_genes

			expr = expr[keep_cols]

			
			gene_list, qv = utils.fetch_geneset(geneset,genesetDir,drug,tissue,fdrThresh,numGenes)
			
			# print(expr.shape)
			
			
			# print(gene_list)
			# print("\n")
			# print(len(gene_list))
			common_genes = [x for x in gene_list if x in expr.columns[1:]]
			# print(len(common_genes))
			# continue
			#for topology in (topProb:=tqdm.tqdm(graphTops, leave=False)):
			for topology in graphTops:
				# topProb.set_description("working on {t} topology".format(t=topology))
		
				networkFile = utils.make_file_path(networkDir,[topology],weighting,".pickle")
		
				resPath = "".join([resultDir,"/".join(["biomarkers","correlation",geneset,drug,tissue, topology]),"/"])
				os.makedirs(resPath,exist_ok = True)
		
				with open(networkFile,"rb") as istream:
					PPI_Graph = pickle.load(istream)
	
			

				issues = list(set(expr.columns[1:]).difference(set(gene_list)))
				# print("Gene list length: {L}".format(L =len(gene_list)))
				# print("Common Genes length: {L}".format(L =len(common_genes)))


				if doOneHop:
					seeds = [x for x in common_genes if x in PPI_Graph.nodes()]
					one_hop = [x for x in seeds]
					for n in seeds:
						nbrs = PPI_Graph.neighbors(n)
						one_hop.extend([x for x in nbrs if x in keep_genes])
					common_genes = [x for x in one_hop]
					# print(common_genes)

			
				LCC_Graph = utils.harmonize_graph_and_geneset(PPI_Graph,common_genes)
				# print("\n")
				# print(LCC_Graph)
				subsetExpression = expr[['Run_ID']+[x for x in LCC_Graph.nodes()]]
			
				data = subsetExpression.merge(resp, on = 'Run_ID')

			
			
				LCC_Graph, gene_to_idx, idx_to_gene = utils.remap_LCC(LCC_Graph,renameField)
				# print(LCC_Graph)
				edgeCurvatures = defaultdict(list)
				nodeCurvatures = defaultdict(list)

				for i in tqdm.tqdm(range(nIters)):
					subsampledData = data.sample(frac = subsamplePct,random_state = rng)
					responderGraph = LCC_Graph.copy()
					nonresponderGraph = LCC_Graph.copy()
					
					responders = subsampledData[subsampledData["Response"]==1].drop(columns = ["Response"])
					nonResponders = subsampledData[subsampledData["Response"]==0].drop(columns = ["Response"])
					
					

					for Graph, dataset, r  in zip([responderGraph, nonresponderGraph], [responders, nonResponders],[1,0]):
					
						genes = nx.get_node_attributes(Graph, renameField)
						total_edge_weight = 0
						for edge in Graph.edges():
							node1, node2 = edge[0], edge[1]
							gene1, gene2 = genes[node1],genes[node2]
							corr,pval =stats.spearmanr(np.log2(dataset[gene1]+1),np.log2(dataset[gene2]+1))
							weight = np.round(0.5*(1+corr),5)
							total_edge_weight += weight
							Graph[node1][node2][weightField] = weight

						if normalizeWeights:
							for edge in Graph.edges():
								node1, node2 = edge[0],edge[1]
								Graph[node1][node2][weightField] = scalingConstant*(Graph[node1][node2][weightField]/total_edge_weight)
							
						nc.assign_densities(Graph,
							alpha = alpha, 
							weight = weightField, 
							measure_name = densityField)
						D = nc.make_APSP_Matrix(Graph,weighted = True, weight=weightField)
						# print("computing curvatures")
						nc.compute_OR_curvature(Graph,D,densityField, edgeCurvatureField, sinkhorn, regParam)
						# print("computing node curvatures")
						nc.compute_node_curvatures(Graph, weight = weightField,
							edge_curv = edgeCurvatureField, 
							node_curv = nodeCurvField,
							norm_node_curv = normNodeCurvField)
					
						for edge in Graph.edges(data=True):
							gene1, gene2 = genes[edge[0]],genes[edge[1]]
							edgeString = gene1+"-"+gene2
							edgeCurvatures['iter'].append(i)
							edgeCurvatures['edge'].append(edgeString)
							edgeCurvatures['ORC'].append(edge[2][edgeCurvatureField])
							edgeCurvatures['Response'].append(r)
					

						for node in Graph.nodes(data=True):
							gene = node[1][renameField]
							nodeCurvatures['iter'].append(i)
							nodeCurvatures['Gene'].append(gene)
							nodeCurvatures['NodeCurve'].append(node[1][nodeCurvField])
							nodeCurvatures['NormNodeCurve'].append(node[1][normNodeCurvField])
							nodeCurvatures['Response'].append(r)
				# print(resPath)	
				statString = "alpha:\t{a}\nfdrThresh:\t{t}\nNumber Unique Genes:\t{g}".format(a=alpha,t=fdrThresh,g=numGenes)
				statString = statString + "\nMinimum TPM:\t{m}".format(m=minTPM)
				statString = statString + "\nSinkhorn:\t{s}".format(s=sinkhorn)
				statString = statString + "\nEpsilon\t{e}".format(e=regParam)
				statString = statString + "\nNum Iters\t{n}".format(n=nIters)
				statString = statString + "\nPct:\t{p}".format(p=subsamplePct)
				statString = statString + "\nQuantile Cutoff:\t{p}".format(p=qv)
				statname = resPath + "stats.txt"
				statname = resPath + "stats.txt"
				with open(statname,"w") as ostream:
					ostream.write(statString)
				nodeDF = pd.DataFrame(nodeCurvatures)
				nodeDF.to_csv(resPath + "node_curvatures.csv")
				edgeDF = pd.DataFrame(edgeCurvatures)
				edgeDF.to_csv(resPath + "edge_curvatures.csv")
				


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)
	
	