# plasma-proteomic-markers-prs-t2d-scripts
 Analysis scripts in support of T2D polygenic score and plasma proteome manuscript 

# README


This repository houses and documents the code used to generate the results in the study Loesch DP *et al.* Identification of plasma proteomic markers underlying polygenic risk of type 2 diabetes and related comorbidities. *medRxiv*, doi: 10.1101/2024.03.15.24304200 (https://www.medrxiv.org/content/10.1101/2024.03.15.24304200v1).

The code found in the analysis_scripts folder has not been designed to regenerate the results as-is. It was written to run on a SLURM based high-performance computing cluster. Job submission scripts were not included, but example job submission scripts can be provided upon request. Note that this project only used publically-availble and previously published software. See the Software section below for a listing of all software. The data analysed in this study is not included in this repository. See the Data Availability section in the preprint for more information.   

# Repository Structure
The scripts archived in analysis_scripts are grouped according to analysis type and are ordered in roughly the same manner as the methods section of the manuscript. 

The sample_pipeline folder contains scripts that can be run on the provided sample data or new data in order to perform an analysis modled after this project. Within this folder, there is a test_data subfolder, a pipeline_scripts subfolder, and a bash script with an example analysis.  

# Software and versions:

Data processing, regression, and survival analyses:
 - R 4.0.02
 - stats R package
 - survival R package (https://github.com/therneau/survival)
Causal inference:
 - R 4.2.2
 - medflex v0.6-10 (https://github.com/jmpsteen/medflex)
 - MendelianRandomization v0.7.0 (https://cran.r-project.org/web/packages/MendelianRandomization/index.html)
 - coloc v5.1.0.1 (https://github.com/chr1swallace/coloc)

Genetic analyses and genotype quality control: 
- PLINK v1.90b6.18 (https://www.cog-genomics.org/plink/2.0/)
- PLINK v2.00a4LM (https://www.cog-genomics.org/plink/2.0/)
- EAGLE 2.4.1 (https://alkesgroup.broadinstitute.org/Eagle/)
- BEAGLE 4.1 (https://faculty.washington.edu/browning/beagle/b4_1.html).
- PRS-CS (https://github.com/getian107/PRScs)
- PRS-CSx (https://github.com/getian107/PRScsx).
- REGENIE v3.1.2 (https://github.com/rgcgithub/regenie).

Pathway analyses was performed using the g:Profiler web portal (https://biit.cs.ut.ee/gprofiler/gost). 

# Sample pipeline
In addition to housing the analysis scripts used for the study, this repository also contains a sample pipeline and a test data set (synthetically generated). This pipeline can be used to perform the analyses used in this study on new data. 
 - Installation: downlaod the sample_pipeline folder to a local machine. The R scripts will install any needed packages. 
 - Usage: all scripts have argument handling with help messages. The sample_pipeline.sh presents the recommended analysis order and additional documentation. The test_data folder contains everything needed to run an example analysis. 
 - Using on new data: note that PGS estimation, QC, and data wrangling needs to occur prior to running the scripts. Please use the sample data, found in test_data, for a model on how the data should be structured.  
