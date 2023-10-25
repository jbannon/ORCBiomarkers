import sys
from typing import List, Dict
import networkx as nx

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
