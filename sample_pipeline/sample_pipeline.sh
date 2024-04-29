#### This script demonstrates a sample pipeline for conducting an analysis like the one found in Loesch et al. 2024. #####
#### Please see the analysis_scripts folder for the actual scripts used the generate the manuscript results ####
# R version 4.3.0 was used for this pipeline 

# specify paths for data and scripts
pipeline=./pipeline_scripts
data=./test_data

#### Part 1: test PGS for association with protein levels #######
# To see if PGS-protein associations are polygenic, the script re-runs the association analysis with pQTL adjustments 
# The script also runs +/- BMI 

Rscript $pipeline/pgs_protein_assoc.R --pheno_file $data/simulated_pheno_data.txt \
--protein_file $data/simulated_protein_data.txt \
--dosage_file $data/simulated_dosage.txt 

# "Rscript $pipeline/pgs_protein_assoc.R --help" to see all options. Script will run using default arguments 


##### Part 2: perform mediation analysis ####
# tests if protein and PGS are association with incident disease/event, and then performs mediation using the medlfex package #

Rscript $pipeline/pgs_protein_mediation.R --pheno_file $data/simulated_pheno_data.txt \
--protein_file $data/simulated_protein_data.txt 


#### Part 3: Cox regression using proteins as exposure #####
# Tests proteins (e.g., PGS-associated proteins ) for assocaiton with time to event, such as time to a clinical trial endpoint or diagnosis
#fits 3 models, one with just standard covariates, another with covariates + clinical risk factors, and a third with covariates + clinical risk factors + biomarkers
Rscript $pipeline/protein_cox_regression.R --pheno_file $data/simulated_pheno_data.txt \
--protein_file $data/simulated_protein_data.txt \
--event_file $data/simulated_time_to_event.txt \
--exclude_file $data/simulated_time_to_event_exclusions.txt \
--biomarkers "BIOMARKER"

##### Part 4: MR and coloc ####### 
#running MR analysis using the MendelianRandomization R package
Rscript $pipeline/run_MR.R --trait_gwas $data/TRAIT_GWAS.txt \
--exposure_gwas $data/PROTEIN_GWAS.txt \
--instruments $data/pqtl_list

#running coloc
Rscript $pipeline/run_coloc.R --trait_gwas $data/TRAIT_GWAS.txt \
--exposure_gwas $data/PROTEIN_GWAS.txt \
--index_snp "rs505922"

### Note that for pathway anaylsis, the gProfiler web portal was used (https://biit.cs.ut.ee/gprofiler/gost) 
# protein/gene lists were uploaded to the server, and the statistical domain was restricted to proteins on the Olink panels 


