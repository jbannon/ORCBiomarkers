import os
import sys
from typing import List






def fetch_genesets(
	data_base_dir:str = "/Users/jamesbannon/Desktop/CancerResponseData/data/processed"
	)->None:
	cmd = "rm ../data/genesets/*"
	os.system(cmd)
	cmd = "cp {f}/genesets/* ../data/genesets".format(f=data_base_dir)
	os.system(cmd)


def fetch_networks(
	data_base_dir:str = "/Users/jamesbannon/Desktop/CancerResponseData/data/processed/networks",
	data_sets:List[str] = ['cri','ccle']
	)->None:
	for dataset in data_sets:
		path = "{p}/{d}".format(p = data_base_dir, d=dataset)
		local_base_path ="../data/networks/{d}".format(d=dataset)

		os.makedirs(local_base_path,exist_ok = True)
		
		cmd = "rm -rf {p}/*".format(p=local_base_path)
		os.system(cmd)
		
		for item in os.listdir(path):
			if item[0]==".":
				continue
			# print("{p}/{i}".format(p=path,i=item))
			item_path = "{p}/{i}".format(p=path,i=item)
			if os.path.isdir(item_path):
				local_item_path = "{p}/{i}".format(p=local_base_path, i = item)
				os.makedirs(local_item_path, exist_ok = True)
				cmd = "cp {p}/{i}/* {pp}/".format(p=path,i=item, pp=local_item_path)
				os.system(cmd)
			else:
				cmd = "cp {p}/{i} {pp}/".format(p=path,i= item, pp= local_base_path)
				os.system(cmd)
		
def fetch_expression(
	data_base_dir:str = "/Users/jamesbannon/Desktop/CancerResponseData/data/processed/expression",
	data_sets:List[str] = ['cri','ccle'],
	file_names:List[str] = ['expression.csv','response.csv','counts.csv']
	)->None:

	for dataset in data_sets:
		path = "{p}/{d}".format(p = data_base_dir, d=dataset)
		local_base_path ="../data/expression/{d}".format(d=dataset)
		os.makedirs(local_base_path,exist_ok = True)
		
		cmd = "rm -rf {p}/*".format(p=local_base_path)
		os.system(cmd)
		

		for item in os.listdir(path):
			
			if item[0]==".":
				continue
			
			item_path = "{p}/{i}".format(p=path,i=item)
			# print(item_path)
			if os.path.isdir(item_path):
				local_item_path = "{p}/{i}".format(p=local_base_path, i = item)
				os.makedirs(local_item_path, exist_ok = True)
				cmd = "rm -rf {p}/*".format(p=local_item_path)
				
				for tissue in os.listdir(item_path):
					if tissue == 'PANCAN':
						continue
					tissue_path = "{p}/{d}/{i}/{t}".format(p=data_base_dir,d=dataset,i=item,t=tissue)
					if os.path.isdir(tissue_path):
						local_tissue_path = "{p}/{i}/{t}/".format(p=local_base_path,i = item,t=tissue)
						os.makedirs(local_tissue_path, exist_ok = True)
						for file in file_names:
							cmd = "cp {p}/{f} {pp}".format(p=tissue_path, f=file,pp = local_tissue_path)
							os.system(cmd)
								

# fetch_genesets()
# fetch_networks()
fetch_expression()