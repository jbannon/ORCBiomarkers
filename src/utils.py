from typing import List
import networkx as nx



def make_file_name(
	components:List[str]
	) -> str:
	# assume components[0] is base_dir, components[-1]
	base_dir = components[0]
	base_dir = base_dir if base_dir[-1]!="/" else base_dir[:-1]
	print(base_dir)

def harmonize_graph_and_geneset(
    G:nx.Graph,
    gene_set:List[str]
    ) -> nx.Graph:
    

    common_genes = list(set(G.nodes).intersection(set(gene_set)))
    # print(len(common_genes))


    # for cc in nx.connected_components(G):
    #     for g in common_genes:
    #         if g in cc:
    #             print("{g}\t{cc}".format(g=g,cc=len(cc)))
    
    H = G.subgraph(common_genes)
    # print(len(H.nodes))
    
    # for cc in nx.connected_components(H):
    #     if len(cc)==1:
    #         print(cc)
    
    if not nx.is_connected(H):
        LCC_genes = max(nx.connected_components(H), key=len)
        H = H.subgraph(LCC_genes)
    
    return H