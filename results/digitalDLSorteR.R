library(digitalDLSorteR)

DDLSObject <- loadRealSCProfiles(
   single.cell.real = sc.chung.breast, ## vector with three elements: single-cell counts, cells metadata and genes metadata
   cell.ID.column = "Cell_ID", gene.ID.column = "external_gene_name",
   min.cells = 1, min.counts = 1, project = "Esa_example" #is it filtering genes with 0 counts?
   )

#--------------------------------------------------------------------------------------------------------------------------
# 
# #Introduce scGAN params?
# library(reticulate)
# 
# #data has to be in the `.h5` or `.h5ad` format
# import scanpy.api as sc
# import os
# import csv
# 
# #Dataset paths
# dataset_dir = "/path_to_where_you_saved/Fresh_68k_PBMCs/hg19/"
# data_file = "matrix.mtx"
# var_names_file = "genes.tsv"
# obs_names_file = "barcodes.tsv"
# output_h5ad_file = "68kPBMCs.h5ad"
# 
# data_path = os.path.join(dataset_dir,data_file)
# var_names_path = os.path.join(dataset_dir,var_names_file)
# obs_names_path = os.path.join(dataset_dir,obs_names_file)
# output_h5ad_path = os.path.join(dataset_dir,output_h5ad_file)
# 
# #Load variable gene names
# with open(var_names_path, "r") as var_file:
#   var_read = csv.reader(var_file, delimiter='\t')
# var_names = []
# for row in var_read:
#   print(row)
# var_names.append(row[1])
# 
# #Load the observations (UMI) names
# with open(obs_names_path, "r") as obs_file:
#   obs_read = csv.reader(obs_file, delimiter='\t')
# obs_names = []
# for row in obs_read:
#   #print(row)
#   obs_names.append(row[0])
# 
# #Load data and transpose it (take a while)
# andata = sc.read(data_path)
# andata = andata.transpose()
# 
# #Load the variable and observation names (and make them unique) into the AnnData
# andata.var_names = var_names
# andata.var_names_make_unique()
# andata.obs_names = obs_names
# andata.obs_names_make_unique()
# 
# #Write the Anndata to hd5ad file
# andata.write(filename=output_h5ad_path)
# 
# 



#Pre-processing (can train sequencially if --train after --process)
#python main.py --param parameters.json --process 

#Training
#python main.py --param parameters.json --train 
#simulate new profiles simSingleCellProfiles

#Generate cells
#python main.py --param parameters.json --generate --cells_no 1000 500 0 200 --model_path path/to/my/model --save_path where_to_save.h5ad


#------------------------------------------------------------------------------------------------------


#Generate composition cell matrix for bulk generation

## data.frame with prior information about cell types
## 3 columns: 
#A cell type column with the same name of the cell type column in cells.metadata. 
#If the name of the column is not the same, function returns an error. Cell types must appear on cells.metadata.
##A second column named 'from' with the start frequency for each cell type.
###A third column named 'to' with the final frequency for each cell type.

 probMatrix <- data.frame(
   Cell_types = c("ER+", "HER2+"),#ER+ will have from 30 to 70
   from = c(rep(30, 4), 1, rep(0, 8)),
   to = c(rep(70, 4), 50, rep(15, 8)) 
   )
 DDLSObject <- generateTrainAndTestBulkProbMatrix(
   object = DDLSObject,
   cell.type.column = "Cell_type",
   prob.design = probMatrix,
   n.bulk.samples = 200 # NUMBER of bulk samples we want to generate. 30000 in real situations
   )
 
 
 ## simulation of bulk samples
  DDLSObject <- generateBulkSamples(
    DDLSObject,
    threads = 2,
    type.data = "both",
    file.backend = "DDLS_bulk.h5" # if NULL, works in-memory
    )
  
  ## preparing samples for training (train) and prediction (test)
  DDLSObject <- prepareDataForTraining(
    object = DDLSObject,
    type.data = "both",
    combine = "bulk", # or both or single-cell
    file.backend = "DDLS_final.h5" # if NULL, works in-memory
    )
  
   ## training using train set and prediction over test set
   DDLSObject <- trainDigitalDLSorterModel(
     object = DDLSObject,
     batch.size = 128,
     num.epochs = 20,
     val = FALSE # if TRUE, a subset of data is used for evaluation during training
     )
   
   
   #METRICS 
  ## calculation of absolute error and squared error
    DDLSObject <- calculateEvalMetrics(DDLSObject)
    ## distribution of absolute errors by cell type
    distErrorPlot(
      DDLSObject,
      error = "AbsErr", # or ppAbsErr, RsqrErr...
      facet.by = "CellType", # or nMix or NULL...
      color.by = "nMix", error.labels = TRUE
    )
    ## correlation between expected and prediction
    corrExpPredPlot(
      DDLSObject, corr = "both", ## pearson or CCC
      facet.by = "CellType", color.by = "CellType",
      filter.sc = FALSE
      )
    ## Bland Altman agreement plot
    blandAltmanLehPlot(
      DDLSObject, facet.by = "CellType",
      color.by = "CellType",
      log.2 = TRUE, density = TRUE
      )
   
#New bulk RNA-Seq deconv
    ## new SummarizedExperiment object
    library(SummarizedExperiment)
     TCGA.breast <- SummarizedExperiment(
       assay = list(counts = TCGA.breast.small)
       )
     ## load SE object into DigitalDLSorter object
     DDLSObject <- loadDeconvDataFromSummarizedExperiment(
       object = DDLSObject,
       se.object = TCGA.breast,
       name.data = "TCGA.breast"
       )
     ## deconvolution using this data
     DDLSObject <- deconvDigitalDLSorterObj(
       object = DDLSObject,
       name.data = "TCGA.breast",
       normalize = TRUE
       )
     ## see results
     barPlotCellTypes(DDLSObject, name.data = "TCGA.breast")
