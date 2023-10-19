import pandas as pd
import numpy as np 
import networkx as nx
import ot







expr = pd.read_csv('../data/expression/Atezo/KIRC/expression.csv')
resp = pd.read_csv('../data/expression/Atezo/KIRC/response.csv')
DE = pd.read_csv("../data/genesets/Atezo_KIRC_DE.csv",index_col = 0)

print(expr)
print(resp)
print(DE)
