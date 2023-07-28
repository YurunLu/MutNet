# MutNet
A heterogeneous network embedding based method by integrating genomic data with genes’ biological processes including protein-protein interactions and pathways.
## Installation
1. Run the following commands for installation:
```
wget https://github.com/AMSSwanglab/MutNet/archive/master.zip  
unzip master.zip
cd MutNet-master
```  
2. Prepare necessary datasets for MutNet including:

Pathway from Reactome:
```
cd ./pathways
wget https://reactome.org/download/current/ReactomePathways.gmt.zip
unzip ReactomePathways.gmt.zip
wget https://reactome.org/download/current/ReactomePathwaysRelation.txt
```
PPI from STRING (and preprocess):
```
cd ./PPI
wget https://stringdb-downloads.org/download/protein.links.v11.5.txt.gz
gunzip protein.links.v11.5.txt.gz
python prePPI.py
```

## Run MutNet
1. Prepare input file in the `./data`: MutNet takes a `.csv` file as input, with which each line shoud have at least 4 columns: `Gene`, `Tumor_Sample`, `Tumor_Type`, `Variant_Classification`. A sample input data can be download from http://karchinlab.org/data/Protocol/pancan-mutation-set-from-Tokheim-2016.txt.gz.
2. Edit `args_profile.txt`, crucial information including:
*  `Mutation_file` for path of input file; 
*  `Pathway_Dir` for path to pathways datasets;
*  `Result_Path` for path to PPI dataset; 
*  `PPI_file` for input file name.
3. After preparing the input files as above, run the following command:
```
python getembed.py
python ConstructNet.py
python outdrivers.py
```
Consturction of network for propagation can be time-consuming. Constructed networks for pan-cancer datasets can be accessed with the following command:
```
cd ./Result/PanCancer/
unrar ./Network_pathwayn.rar
unrar ./Network_ppi.rar
```
and then `python ConstructNet.py` can be skipped.

4. A fold will be constructed in the `Result_Path` fold according to the args, and predicted cancer gene will be listed in `driver_genes.csv` in the fold.

A list for all genes with MutScore_raw and MutScore_AfterNP is also provided in `gene_info_proped.csv`

## Requirements
Python environment: python 3

Python packages: 
*  adjustText==0.8
*  gensim==3.8.3
*  matplotlib==3.5.3
*  matplotlib_venn==0.11.6
*  networkx==2.6.3
*  numpy==1.21.6
*  pandas==1.3.5
*  preprocess==2.0.0
*  scikit_learn==1.0.2
*  scipy==1.7.3
*  seaborn==0.12.2
*  matplotlib==3.3.2

Memory: >= 3.0 Gb

## Citation

If you use MutNet or MutNet associated resources, please cite

Yurun Lu, et al. Integrative analysis of pan cancer genomic data with heterogeneous network representation learning. 2023.

Codes for reproducing results in the article are also included in `Codes for Other Results`.


