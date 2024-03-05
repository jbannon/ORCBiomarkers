import pandas as pd 



df = pd.read_csv("../data/expression/cri/Ipi+Pembro/SKCM/response.csv")
print(df)
print(len(pd.unique(df['Run_ID'])))