# UTRExtension
# Identification Untranslated Region (UTR) Extensions/Variations in RNA-seq data

## Identification of UTR Extensions
* A computational pipeline was designed to identify 5’/3’ UTR extensions in RNA-seq data (Figure 1).

![Schematic representation of computational pipeline for identification 5’/3' UTR extensions.](https://github.rcac.purdue.edu/BioinformaticsCore/UTRExtension/blob/master/workflow.png)
* Using 10 bp window and using 5’ start and 3’ end of each gene the 5’/3’ extensions were contiguously assembled on either side of the gene. 
* The window of 10bp was extended to next 10 bp if the RPKM >=1. 
* The RPKM values of these windows were calculated if the reads mapped exclusively to the intergenic regions, the reads partially mapped to the genic and intergenic regions were not considered in the RPKM calculation. 
* Additionally, only the uniquely mapped reads with mapping quality greater than 50 were considered. Using this technique the 5’/3’ UTR extensions were predicted from all the samples and later all extensions were combined together while keeping the longest extension per gene. 
* After obtaining the extensions, the expression of the assembled extensions was determined under stress and control conditions.
* The final extensions were retained if the absolute log 2 fold change was greater than 1 in stress condition than normal using the DESeq2  method at FDR <=0.05. 


## Run it from command line

module load anaconda/5.1.0-py27  
conda env create -f environment.yaml  
## This will create the conda environment for the pipeline
conda activate snakemake_py2  
## Command for the pipeline
snakemake -s utr_pipeline.v1.py --configfile config.yaml -j 30 --cluster-config cluster.json --cluster "qsub -q {cluster.queue} -l nodes=1:ppn={cluster.n} -l walltime={cluster.time} -e {cluster.error} -o {cluster.output} -N {cluster.name} -M {cluster.email}"  > utr_pipeline.log  2>&1 &

Paper followed: https://genome.cshlp.org/content/27/8/1427.long 
