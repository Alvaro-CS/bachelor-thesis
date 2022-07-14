import scanpy as sc
import os
import csv

dataset_dir = ""
data_file = "counts.mtx"
var_names_file = "features.tsv"
obs_names_file = "barcodes.tsv"
output_h5ad_file = "pesaOUT.h5ad"

data_path = os.path.join(dataset_dir,data_file)
var_names_path = os.path.join(dataset_dir,var_names_file)
obs_names_path = os.path.join(dataset_dir,obs_names_file)
output_h5ad_path = os.path.join(dataset_dir,output_h5ad_file)

with open(var_names_path, "r") as var_file:
    var_read = csv.reader(var_file, delimiter='\t')
    var_names = []
    for row in var_read:
        #print(row)
        var_names.append(row[1])

with open(obs_names_path, "r") as obs_file:
    obs_read = csv.reader(obs_file, delimiter='\t')
    obs_names = []
    for row in obs_read:
        #print(row)
        obs_names.append(row[0])

andata = sc.read(data_path) 
andata = andata.transpose()

andata.var_names = var_names
andata.var_names_make_unique()
andata.obs_names = obs_names
andata.obs_names_make_unique()

andata.write(filename=output_h5ad_path)