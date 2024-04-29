###wrapper script for running reverse MR #####
protein=$1

export PATH="../envs/MR/bin/:$PATH"

#check if input file exists
#if [ ! -f ./data/reverse_mr_data/inputs/NASH_AST_ALT_${protein}_reverse_MR_inputs.tsv ]; then
    echo "preparing inputs"
    Rscript ./scripts/06_coloc_MR/reverse_MR/prep_reverse_MR_data.R $protein
    echo "done preparing inputs"
#fi 

#run reverse MR
Rscript ./scripts/06_coloc_MR/reverse_MR/run_reverse_MR.R $protein