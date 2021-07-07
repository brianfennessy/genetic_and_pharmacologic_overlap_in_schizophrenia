# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Store user home directory and navigate to it
user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
setwd(user_home_directory)

# Get the official home directory of the user
home_dir <- getwd()

# Navigate to the raw data
setwd(home_dir)

source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")
tar[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load APs and PPs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/propsychotics_final.tsv")
ap <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/ap_drugs_with_targets.tsv")

pp[,rxcui := as.character(rxcui)]
ap[,rxcui := as.character(rxcui)]

pp[! rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name
ap[! rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Isolate targets of each respective drug set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_targs <- unique(tar[rxcui %in% pp$rxcui,]$gene)
ap_targs <- unique(tar[rxcui %in% ap$rxcui,]$gene)

length(pp_targs) # [1] 963
length(ap_targs) # [1] 433
length(intersect(pp_targs,ap_targs)) # [1] 312

pp_targs_final <- as.data.table(pp_targs)
ap_targs_final <- as.data.table(ap_targs)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(pp_targs_final, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/propsychotic_targets_all.tsv")
fwrite(ap_targs_final, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/antipsychotic_targets_all.tsv")