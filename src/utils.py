import numpy as np
from scipy.stats import kurtosis,mode, skew, beta
import sys
from typing import List, Dict, Union
import networkx as nx
import pandas as pd
from collections import defaultdict

DRUG_TISSUE_MAP = {"Atezo":["KIRC","BLCA"],"Pembro":["SKCM","STAD"],
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
    pathNames:List[str], 
    fname:str, 
    ext:str
    )->str:
    
    pathNames.append("")
    ext = ext if ext[0]=="." else "." + ext

    base_dir = base_dir[:-1] if base_dir[-1]=="/" else base_dir
    path = [base_dir]
    path.extend(pathNames)
    path ="/".join([x for x in path])
    
    file_path = "".join([path,fname,ext])
 
    

    return file_path

def harmonize_graph_and_geneset(
    G:nx.Graph,
    gene_set:List[str]
    ) -> nx.Graph:
    

    common_genes = [x for x in list(G.nodes) if x in gene_set]



    H = G.subgraph(common_genes)
    
    if len(H.nodes)>0:
        if not nx.is_connected(H):
            
            LCC_genes = max(nx.connected_components(H), key=len)
            H = H.subgraph(LCC_genes)
    
    return H


def remap_LCC(
    LCC_Graph:nx.Graph,
    newFieldName:str = 'Gene'
    ):
    gene_to_idx  = {} 
    idx_to_gene = {}
    for idx, gene in enumerate(LCC_Graph.nodes):
        gene_to_idx[gene] = idx
        idx_to_gene[idx] = gene
    
    LCC_Graph = nx.relabel_nodes(LCC_Graph,gene_to_idx)
    nx.set_node_attributes(LCC_Graph,idx_to_gene, newFieldName)
    return LCC_Graph, gene_to_idx, idx_to_gene




def empiricalBayesGeneSelection(
    df:pd.DataFrame,
    probThresh:float,
    cutType:str,
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

    



