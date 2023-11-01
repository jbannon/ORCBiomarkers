import sys
from typing import List, Dict
import networkx as nx
import pandas as pd


DRUG_TISSUE_MAP = {"Atezo":["KIRC","BLCA"],"Pembro":["SKCM","STAD"],
    "Nivo":["KIRC","SKCM"], "Ipi":["SKCM"], "Ipi+Pembro":["SKCM"]}



def fetch_geneset(
    geneset:str,
    gsDir:str,
    drug:str, 
    tissue:str,
    fdrThresh:float,
    nGenes:int = 150
    ) -> List[str]:
    if geneset.upper()=="DE":
        diffExpFile = make_file_path(gsDir,[],"{d}_{t}_DE".format(d=drug, t=tissue),".csv")
        DE = pd.read_csv(diffExpFile)
        DE = DE[DE['Thresh.Value']==fdrThresh]
        
        q = (1.0*nGenes)/len(pd.unique(DE['Gene']))
        thresh = 1-q
        
        qval = DE['Count'].quantile(thresh)
        DE =  DE[DE['Count'] >= qval]
        gene_list = list(DE['Gene'].values)
        
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
    # print(len(common_genes))


    H = G.subgraph(common_genes)
    

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