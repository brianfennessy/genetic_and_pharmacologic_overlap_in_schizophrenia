# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET - libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Store user home directory and navigate to it
user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
setwd(user_home_directory)

# Get the official home directory of the user
home_dir <- getwd()

# Navigate to the raw data
setwd(home_dir)

source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - MedDRA psychosis terms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meddra <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/meddra_psychosis.tsv")
term_names <- unique(meddra[,.(pt_code,pt_name)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")
tar[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - RxNorm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rxn <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/rxnorm_refined.tsv")
rxn_in <- unique(rxn[tty=="IN",])
rxn_names <- unique(rxn_in[sab=="RXNORM",.(rxcui,drug_name)])
rxn_names[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - ATC to RxNorm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atc_rxn <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/atc_rxnorm_mapped.tsv")
atc_desc <- unique(atc_rxn[,.(rxcui,level3_description)])
atc_desc[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - Antipsychotics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ap <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/ap_drugs.tsv")
ap[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - SIDER
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdr_master <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/sider_original.tsv")
uniqueN(sdr_master$rxcui) # [1] 1332

uniqueN(sdr_master[,.(rxcui,meddra)]) # [1] 147959
uniqueN(sdr_master$rxcui) # [1] 1332
uniqueN(sdr_master$meddra) # [1] 4806

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LIMIT - SIDER to only adverse event rows and adverse events attributed to psychosis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gold_sdr <- unique(sdr_master[meddra_types=="ae",])
gold_sdr <- gold_sdr[,.(rxcui,pt_code=meddra,meddra_level)]

uniqueN(gold_sdr[,.(rxcui,pt_code)]) # [1] 133888
uniqueN(gold_sdr$rxcui) # [1] 1298
uniqueN(gold_sdr$pt_code) # [1] 4201

sdr_pp <- unique(gold_sdr[pt_code %in% meddra$pt_code,])
sdr_pp[,rxcui:=as.character(rxcui)]

uniqueN(sdr_pp[,.(rxcui,pt_code)]) # [1] 1143
uniqueN(sdr_pp$rxcui) # [1] 433
uniqueN(sdr_pp$pt_code) # [1] 46

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE - Any rxcui's in our SIDER prospychotics that is found in our list of antipsychotics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdr_pp <- sdr_pp[!(rxcui %in% ap$rxcui),.(rxcui,pt_code)]

uniqueN(sdr_pp[,.(rxcui,pt_code)]) # [1] 1080
uniqueN(sdr_pp$rxcui) # [1] 412
uniqueN(sdr_pp$pt_code) # [1] 43

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ADD - ATC, MedDRA, and RxNorm information to SIDER propsychotics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdr_pp <- merge(sdr_pp,term_names,by="pt_code",all.x=TRUE)
sdr_pp <- merge(sdr_pp,rxn_names,by="rxcui",all.x=TRUE)
sdr_pp <- merge(sdr_pp,atc_desc,by="rxcui",all.x=TRUE)
sdr_pp[is.na(level3_description),level3_description := "UNMAPPED"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ISOLATE - Sets of SIDER drugs and SIDER propsychotics that have target data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdr_pp_w_targs <- sdr_pp[rxcui %in% tar$rxcui,]
sdr_master_w_targs <- sdr_master[rxcui %in% tar$rxcui,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE - Sets of SIDER drugs and SIDER propsychotics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(sdr_master_w_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/sider_original_with_targets.tsv")
fwrite(sdr_pp, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/sider_propsychotics.tsv")
fwrite(sdr_pp_w_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/sider_propsychotics_with_targets.tsv")
