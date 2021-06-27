# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Store user home directory and navigate to it
user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
setwd(user_home_directory)

# Get the official home directory of the user
home_dir <- getwd()

# Navigate to the raw data zip file and unzip
setwd(home_dir)

source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- unique(fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")[,.(rxcui,drug_name,gene)])
tar[,rxcui := as.character(rxcui)]
genes <- unique(tar$gene)
num_drugs_in_targets <- uniqueN(tar$rxcui)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load APs and PPs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_all <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_propsychotics_drugs_only.tsv")
pp_int <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/propsychotics_final.tsv")
ap <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/ap_drugs.tsv")
pp_no_cns <- unique(fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/atc_rxnorm_mapped.tsv")[,.(rxcui,drug_name,level1_description)])
pp_no_cns <- pp_no_cns[rxcui %in% pp_all$rxcui & level1_description != "NERVOUS SYSTEM",]
pp_no_cns <- unique(pp_no_cns[,.(rxcui,drug_name)])

pp_all[,rxcui := as.character(rxcui)]
pp_int[,rxcui := as.character(rxcui)]
ap[,rxcui := as.character(rxcui)]
pp_no_cns[,rxcui := as.character(rxcui)]

pp_all <- pp_all[rxcui %in% tar$rxcui,]
pp_int <- pp_int[rxcui %in% tar$rxcui,]
ap <- ap[rxcui %in% tar$rxcui,]
pp_no_cns <- pp_no_cns[rxcui %in% tar$rxcui,]

pp_all[! rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name
pp_int[! rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name
ap[! rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name
pp_no_cns[! rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name

num_pp_all_observed <- uniqueN(pp_all$rxcui)
num_pp_int_observed <- uniqueN(pp_int$rxcui)
num_ap_observed <- uniqueN(ap$rxcui)
num_pp_no_cns_observed <- uniqueN(pp_no_cns$rxcui)

num_pp_all_observed # [1] 240
num_pp_int_observed # [1] 147
num_ap_observed # [1] 46
num_pp_no_cns_observed # [1] 127

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load PP random sampling (all)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load
vgb <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_filter_with_targets_drugs_only.tsv")
vgb[,rxcui := as.character(rxcui)]

# Make sure all drugs have targets in DrugBank or SeaChange
vgb[!rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name

# Remove all antipsychotic drugs
vgb <- vgb[!rxcui %in% ap$rxcui,]

# Save and count
pp_all_random_sampling <- unique(vgb$rxcui)
length(pp_all_random_sampling) # [1] 1951

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load PP random sampling (intersect)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load
sdr <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/sider_original_with_targets.tsv")
sdr[,rxcui := as.character(rxcui)]

# Make sure all drugs have targets in DrugBank or SeaChange
sdr[!rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 4 cols): rxcui,meddra,meddra_level,meddra_types

# Remove all antipsychotic drugs
sdr <- sdr[!rxcui %in% ap$rxcui,]

# Save and count
sdr_drugs <- unique(sdr$rxcui)
length(sdr_drugs) # [1] 1123

# Find intersecting drugs
pp_int_random_sampling <- unique(intersect(pp_all_random_sampling,sdr_drugs))
length(pp_int_random_sampling) # [1] 1106

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load AP random sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load
atc <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/atc_rxnorm_mapped.tsv")

# Make sure all drugs have targets in DrugBank or SeaChange
atc_w_targs <- atc[rxcui %in% tar$rxcui,]

# Remove all propsychotic drugs
atc_w_targs <- atc_w_targs[!rxcui %in% pp_all$rxcui & !rxcui %in% pp_int$rxcui,]

# Save and count
ap_random_sampling <- unique(atc_w_targs$rxcui)
length(ap_random_sampling) # [1] 1403

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load PP random sampling (no CNS)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load
vgb_no_cns <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_filter_with_targets_drugs_only.tsv")
vgb_no_cns[,rxcui := as.character(rxcui)]
cns_drugs <- unique(fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/atc_rxnorm_mapped.tsv")[,.(rxcui,drug_name,level1_description)])
cns_drugs <- cns_drugs[level1_description == "NERVOUS SYSTEM",]

# Make sure all drugs have targets in DrugBank or SeaChange
vgb_no_cns[!rxcui %in% tar$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name

# Remove all antipsychotic drugs and CNS drugs
vgb_no_cns <- vgb_no_cns[!rxcui %in% ap$rxcui,]
vgb_no_cns <- vgb_no_cns[!rxcui %in% cns_drugs$rxcui,]

# Save and count
pp_no_cns_random_sampling <- unique(vgb_no_cns$rxcui)
length(pp_no_cns_random_sampling) # [1] 1697

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Isolate targets of each respective drug set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_targs_all <- unique(tar[rxcui %in% pp_all$rxcui,]$gene)
pp_targs_int <- unique(tar[rxcui %in% pp_int$rxcui,]$gene)
ap_targs <- unique(tar[rxcui %in% ap$rxcui,]$gene)
pp_targs_no_cns <- unique(tar[rxcui %in% pp_no_cns$rxcui,]$gene)

length(pp_targs_all) # [1] 1134
length(pp_targs_int) # [1] 963
length(ap_targs) # [1] 433
length(pp_targs_no_cns) # [1] 916
length(intersect(pp_targs_all,ap_targs)) # [1] 334
length(intersect(pp_targs_int,ap_targs)) # [1] 312
length(intersect(pp_targs_no_cns,ap_targs)) # [1] 296

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set number of permutations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perms <- 100000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random sampling - AP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up the matrix
mymtx_ap <- matrix(0, nrow=length(genes), ncol=(perms+4))
rownames(mymtx_ap) <- genes

# Add observed proportions to the first column
cur <- tar[rxcui %in% unique(ap$rxcui)][,list(frac=.N/num_ap_observed),gene]
mymtx_ap[cur$gene,1] <- cur$frac

# Add random sample proportions
for (i in 1:perms) {
  ap_random_set <- sample(ap_random_sampling, num_ap_observed, replace=FALSE, prob=NULL)
  cur <- tar[rxcui %in% ap_random_set][,list(frac=.N/num_ap_observed),gene]
  mymtx_ap[cur$gene,(i+1)] <- cur$frac
}

# Add mean random proportion for each gene
mymtx_ap[,(perms+2)] <- apply(mymtx_ap[,2:(perms+1)], 1, mean)

# Calculate p-value for each gene
mymtx_ap[,(perms+3)] <- apply(mymtx_ap, 1, function (x) length(which(x[2:(perms+1)] >= x[1]))/(perms))

# Add total proportions to the last column
cur <- tar[,list(frac=.N/num_drugs_in_targets),gene]
mymtx_ap[cur$gene,(perms+4)] <- cur$frac

# Check
mymtx_ap[1:5,c(1,(perms+2),(perms+3),(perms+4))]

# Save
save(mymtx_ap, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/ap_significant_targets.Rdata")

# Check that it loads
# load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/ap_significant_targets.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random sampling - PP (all)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up the matrix
mymtx_pp_all <- matrix(0, nrow=length(genes), ncol=(perms+4))
rownames(mymtx_pp_all) <- genes

# Add observed proportions to the first column
cur <- tar[rxcui %in% unique(pp_all$rxcui)][,list(frac=.N/num_pp_all_observed),gene]
mymtx_pp_all[cur$gene,1] <- cur$frac

# Add random sample proportions
for (i in 1:perms) {
  pp_all_random_set <- sample(pp_all_random_sampling, num_pp_all_observed, replace=FALSE, prob=NULL)
  cur <- tar[rxcui %in% pp_all_random_set][,list(frac=.N/num_pp_all_observed),gene]
  mymtx_pp_all[cur$gene,(i+1)] <- cur$frac
}

# Add mean random proportion for each gene
mymtx_pp_all[,(perms+2)] <- apply(mymtx_pp_all[,2:(perms+1)], 1, mean)

# Calculate p-value for each gene
mymtx_pp_all[,(perms+3)] <- apply(mymtx_pp_all, 1, function (x) length(which(x[2:(perms+1)] >= x[1]))/(perms))

# Add total proportions to the last column
cur <- tar[,list(frac=.N/num_drugs_in_targets),gene]
mymtx_pp_all[cur$gene,(perms+4)] <- cur$frac

# Check
mymtx_pp_all[1:5,c(1,(perms+2),(perms+3),(perms+4))]

# Save
save(mymtx_pp_all, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_all_significant_targets.Rdata")

# Check that it loads
# load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_all_significant_targets.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random sampling - PP (int)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up the matrix
mymtx_pp_int <- matrix(0, nrow=length(genes), ncol=(perms+4))
rownames(mymtx_pp_int) <- genes

# Add observed proportions to the first column
cur <- tar[rxcui %in% unique(pp_int$rxcui)][,list(frac=.N/num_pp_int_observed),gene]
mymtx_pp_int[cur$gene,1] <- cur$frac

# Add random sample proportions
for (i in 1:perms) {
  pp_int_random_set <- sample(pp_int_random_sampling, num_pp_int_observed, replace=FALSE, prob=NULL)
  cur <- tar[rxcui %in% pp_int_random_set][,list(frac=.N/num_pp_int_observed),gene]
  mymtx_pp_int[cur$gene,(i+1)] <- cur$frac
}

# Add mean random proportion for each gene
mymtx_pp_int[,(perms+2)] <- apply(mymtx_pp_int[,2:(perms+1)], 1, mean)

# Calculate p-value for each gene
mymtx_pp_int[,(perms+3)] <- apply(mymtx_pp_int, 1, function (x) length(which(x[2:(perms+1)] >= x[1]))/(perms))

# Add total proportions to the last column
cur <- tar[,list(frac=.N/num_drugs_in_targets),gene]
mymtx_pp_int[cur$gene,(perms+4)] <- cur$frac

# Check
mymtx_pp_int[1:5,c(1,(perms+2),(perms+3),(perms+4))]

# Save
save(mymtx_pp_int, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_int_significant_targets.Rdata")

# Check that it loads
# load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_int_significant_targets.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random sampling - PP (no CNS)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up the matrix
mymtx_pp_no_cns <- matrix(0, nrow=length(genes), ncol=(perms+4))
rownames(mymtx_pp_no_cns) <- genes

# Add observed proportions to the first column
cur <- tar[rxcui %in% unique(pp_no_cns$rxcui)][,list(frac=.N/num_pp_no_cns_observed),gene]
mymtx_pp_no_cns[cur$gene,1] <- cur$frac

# Add random sample proportions
for (i in 1:perms) {
  pp_no_cns_random_set <- sample(pp_no_cns_random_sampling, num_pp_no_cns_observed, replace=FALSE, prob=NULL)
  cur <- tar[rxcui %in% pp_no_cns_random_set][,list(frac=.N/num_pp_no_cns_observed),gene]
  mymtx_pp_no_cns[cur$gene,(i+1)] <- cur$frac
}

# Add mean random proportion for each gene
mymtx_pp_no_cns[,(perms+2)] <- apply(mymtx_pp_no_cns[,2:(perms+1)], 1, mean)

# Calculate p-value for each gene
mymtx_pp_no_cns[,(perms+3)] <- apply(mymtx_pp_no_cns, 1, function (x) length(which(x[2:(perms+1)] >= x[1]))/(perms))

# Add total proportions to the last column
cur <- tar[,list(frac=.N/num_drugs_in_targets),gene]
mymtx_pp_no_cns[cur$gene,(perms+4)] <- cur$frac

# Check
mymtx_pp_no_cns[1:5,c(1,(perms+2),(perms+3),(perms+4))]

# Save
save(mymtx_pp_no_cns, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_no_cns_significant_targets.Rdata")

# Check that it loads
# load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_no_cns_significant_targets.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract relevant data and convert to data tables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pp_all_sig_targs <- data.table(gene=rownames(mymtx_pp_all),
								total_prop=mymtx_pp_all[,(perms+4)],
								observed_prop=mymtx_pp_all[,1],
								mean_rand_prop=mymtx_pp_all[,(perms+2)],
								p_value=mymtx_pp_all[,(perms+3)])

pp_int_sig_targs <- data.table(gene=rownames(mymtx_pp_int),
								total_prop=mymtx_pp_int[,(perms+4)],
								observed_prop=mymtx_pp_int[,1],
								mean_rand_prop=mymtx_pp_int[,(perms+2)],
								p_value=mymtx_pp_int[,(perms+3)])

ap_sig_targs <- data.table(gene=rownames(mymtx_ap),
								total_prop=mymtx_ap[,(perms+4)],
								observed_prop=mymtx_ap[,1],
								mean_rand_prop=mymtx_ap[,(perms+2)],
								p_value=mymtx_ap[,(perms+3)])

pp_no_cns_sig_targs <- data.table(gene=rownames(mymtx_pp_no_cns),
								total_prop=mymtx_pp_no_cns[,(perms+4)],
								observed_prop=mymtx_pp_no_cns[,1],
								mean_rand_prop=mymtx_pp_no_cns[,(perms+2)],
								p_value=mymtx_pp_no_cns[,(perms+3)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bonferroni correction
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Establish corrected p-value
bonferroni_corrected_p <- (0.05/uniqueN(tar$gene))
bonferroni_corrected_p # [1] 1.469724e-05

# Save
pp_all_sig_targs_bonferroni_correction <- pp_all_sig_targs[p_value < bonferroni_corrected_p,]
pp_int_sig_targs_bonferroni_correction <- pp_int_sig_targs[p_value < bonferroni_corrected_p,]
ap_sig_targs_bonferroni_correction <- ap_sig_targs[p_value < bonferroni_corrected_p,]
pp_no_cns_sig_targs_bonferroni_correction <- pp_no_cns_sig_targs[p_value < bonferroni_corrected_p,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Isolate targets with p-value < 0.05
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_all_sig_targs_p_05 <- pp_all_sig_targs[p_value < 0.05,]
pp_int_sig_targs_p_05 <- pp_int_sig_targs[p_value < 0.05,]
ap_sig_targs_p_05 <- ap_sig_targs[p_value < 0.05,]
pp_no_cns_sig_targs_p_05 <- pp_no_cns_sig_targs[p_value < 0.05,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check counts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check counts at p-value less than 0.05
uniqueN(pp_all_sig_targs_p_05$gene) # [1] 169
uniqueN(pp_int_sig_targs_p_05$gene) # [1] 118
uniqueN(ap_sig_targs_p_05$gene) # [1] 129
uniqueN(pp_no_cns_sig_targs_p_05$gene) # [1] 118
uniqueN(intersect(pp_all_sig_targs_p_05$gene,ap_sig_targs_p_05$gene)) # [1] 67
uniqueN(intersect(pp_int_sig_targs_p_05$gene,ap_sig_targs_p_05$gene)) # [1] 52
uniqueN(intersect(pp_no_cns_sig_targs_p_05$gene,ap_sig_targs_p_05$gene)) # [1] 51
 
# Check counts at corrected p-value
uniqueN(pp_all_sig_targs[p_value < bonferroni_corrected_p,]) # [1] 49
uniqueN(pp_int_sig_targs[p_value < bonferroni_corrected_p,]) # [1] 8
uniqueN(ap_sig_targs[p_value < bonferroni_corrected_p,]) # [1] 58
uniqueN(pp_no_cns_sig_targs[p_value < bonferroni_corrected_p,]) # [1] 17

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(pp_all_sig_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_all.tsv")
fwrite(pp_int_sig_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_int.tsv")
fwrite(ap_sig_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/antipsychotic_significant_targets.tsv")
fwrite(pp_no_cns_sig_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_no_cns.tsv")


