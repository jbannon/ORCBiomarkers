import numpy as np
from scipy.stats import kurtosis,mode, skew, beta
import sys
from typing import List, Dict, Union
import networkx as nx
import pandas as pd
from collections import defaultdict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC


DRUG_TISSUE_MAP = {"Atezo":["KIRC","BLCA"],"Pembro":["STAD"],
    "Nivo":["KIRC","SKCM"], "Ipi":["SKCM"], "Ipi+Pembro":["SKCM"],
    "erlotinib":['LUAD'],"crizotinib":['LUAD'],'sorafenib':["LUAD"],'sunitinib':["LUAD"]}


DRUG_DATASET_MAP = {
    'sorafenib':'ccle',
    'erlotinib':'ccle',
    'crizotinib':'ccle',
    'sunitinib':'ccle',
    'Nivo':'cri',
    'Ipi':'cri',
    'Pembro':'cri',
    'Atezo':'cri',
    'Ipi+Pembro':'cri'
}


DRUG_TARGET_MAP = {'Atezo':'PD-L1','Pembro':'PD1','Nivo':'PD1','Ipi':'CTLA4'}


TARGET_GENE_MAP = {'PD-L1':'CD274', 'PD1':'PDCD1', 'CTLA4':'CTLA4'}




def fetch_drug_targets(
    drug:str
    ) -> List[str]:
    if drug in DRUG_TARGET_MAP.keys():
        targets = [TARGET_GENE_MAP[DRUG_TARGET_MAP[drug]]]
    elif drug == "Ipi+Pembro":
        t1 = TARGET_GENE_MAP[DRUG_TARGET_MAP["Ipi"]]
        t2 = TARGET_GENE_MAP[DRUG_TARGET_MAP["Pembro"]]
        targets = [t1,t2]
    else:
        fname = "../data/genesets/{d}_targets.txt".format(d=drug)
        with open(fname, "r") as istream:
            lines = istream.readlines()
        targets = [x.rstrip() for x in lines]

    return targets


def process_pvalue_data(
    pval_df:pd.DataFrame,
    gene_col:str, 
    pval_thresh:float,
    pval_col:str = 'adj_pvals',
    )->List[str]:
    

    
    pval_df = pval_df[pval_df[pval_col]<=pval_thresh]

    
    genes = pval_df[gene_col].values
    gene_list = []
    

    for gene in genes:
        nodes = [x for x in gene.split(";")]
        # we split if it's an edge
        for g in nodes:
            gene_list.append(g)

    gene_list = list(pd.unique(gene_list))
    return gene_list
    
def fetch_geneset(
    geneset:str,
    gsDir:str,
    drug:str, 
    tissue:str,
    fdrThresh:float,
    cutoff:Union[int,str] = 'EB',
    probThresh:float = 0.66
    ) -> List[str]:
    


    if isinstance(cutoff,str):
        assert cutoff.upper() in ['EB','EBCDF','QCUT'], "cutoff.upper() must be in ['EB','EBCDF','QCUT']"
    
    elif isinstance(cutoff,int):
        assert cutoff > 0, "cutoff value must be non-negative"

    if geneset.upper()=="DE":
       
        diffExpFile = make_file_path(gsDir,[],"{d}_{t}_DE".format(d=drug, t=tissue),".csv")
        DE = pd.read_csv(diffExpFile)
        DE = DE[DE['Thresh.Value']==fdrThresh]
        
        if isinstance(cutoff,str):
            if cutoff.upper() == 'QCUT':
                qval = DE['Count'].quantile(probThresh)
                DE =  DE[DE['Count'] >= qval]
                gene_list = list(DE['Gene'].values)
                qval = None
            else:
                gene_list = empiricalBayesGeneSelection(DE,probThresh,cutoff)
                qval = None


        else:
            DE = DE[DE['Count']>=cutoff]
            gene_list = list(DE['Gene'].values)
            qval  = None
    else:
        fname = gsDir+geneset.upper()+".txt"
        with open(fname, 'r') as istream:
            lines = istream.readlines()
            lines = [x.rstrip() for x in lines]
        gene_list = lines
        qval = None
    
    return gene_list, qval



def unpack_parameters(
    D:Dict
    ):
    if len(D.values())>1:
        return tuple(D.values())
    else:
        return tuple(D.values())[0]

def make_file_path(
    base_dir:str, 
    path_names:List[str], 
    fname:str, 
    ext:str
    )->str:
    
    path_names.append("")
    ext = ext if ext[0]=="." else "." + ext

    base_dir = base_dir[:-1] if base_dir[-1]=="/" else base_dir
    path = [base_dir]
    path.extend(path_names)
    path ="/".join([x for x in path])
    
    file_path = "".join([path,fname,ext])
 
    

    return file_path

def harmonize_graph_and_geneset(
    G:nx.Graph,
    gene_set:List[str]
    ) -> nx.Graph:
    

    
    common_genes = [x for x in list(G.nodes) if x in gene_set]
    
    G.remove_nodes_from([n for n in G.nodes if n not in common_genes])

    LCC_genes = sorted(list(max(nx.connected_components(G), key=len)))
    G.remove_nodes_from([n for n in G.nodes if n not in LCC_genes])
    
    
    
    return G


def rename_nodes(
    G:nx.Graph,
    new_field_name:str = 'Gene'
    ):
    gene_to_idx  = {} 
    idx_to_gene = {}
    for idx, gene in enumerate(G.nodes):
        gene_to_idx[gene] = idx
        idx_to_gene[idx] = gene
    G = nx.relabel_nodes(G,gene_to_idx)
    nx.set_node_attributes(G,idx_to_gene, new_field_name)
    return G, gene_to_idx, idx_to_gene




def empirical_bayes_gene_selection(
    df:pd.DataFrame,
    probThresh:float,
    cutType:str = "EB",
    conf:float = 0.9,
    nRuns:int = 200):
    

    temp = df.copy(deep=True)
    
    temp["Hits"] = temp['Count']
    temp['Misses'] = nRuns-temp['Count']
    temp['Proportion'] = temp["Count"]/nRuns

    a, b,loc, scale = beta.fit(temp[temp["Proportion"]<1.0]["Proportion"].values,fscale = 1, floc = 0)

    temp["EB_avg"] = (temp["Hits"]+a)/(nRuns+a+b)
            
    EB_point = np.amin(temp[temp["EB_avg"]>=probThresh]["Count"].values)
            
            
    prob_estimates = defaultdict(list)
    # map CDFs
    for idx, row in temp.iterrows():            
        a1 = a + row["Hits"]
        b1 = b + row["Misses"]
        thisBeta = beta(a=a1,b=b1)
        thisProb = thisBeta.sf(probThresh)
        prob_estimates['Gene'].append(row['Gene'])
        prob_estimates['EB_Prob'].append(thisProb)
    pb = pd.DataFrame(prob_estimates)
    EB = pb[pb['EB_Prob']>=conf]
    EB = temp[temp["Gene"].isin(EB["Gene"].values)]

    if cutType == 'EB':
        temp = temp[temp['Count'] >= EB_point]
    elif cutType == "EBCDF":
        try:
            cut = np.amin(EB["Count"].values)
            temp = temp[temp['Count'] >= cut]
        except: 
            print("cdf causing an issue, using prob estimate")
            temp = temp[temp['Count'] >= EB_point]

    gene_list = list(temp['Gene'].values)
    return gene_list

    



def make_model_and_param_grid(
    model_name:str,
    reg_min:float,
    reg_max:float,
    reg_step:float,
    model_max_iters:int
    ):
    
    preproc = ('preproc',StandardScaler())
    
    if model_name == 'LogisticRegression':
        classifier = ('clf',LogisticRegression(class_weight = 'balanced',max_iter = model_max_iters))
    elif model_name == 'LinearSVM':
        classifier = ('clf',LinearSVC(class_weight = 'balanced',max_iter = model_max_iters))
    

    param_grid = {'clf__C':np.arange(reg_min,reg_max,reg_step)}                            
    
    model = Pipeline([preproc,classifier])

    return model, param_grid

