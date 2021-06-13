# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET - Libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - PP and AP significant targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_all_sig_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/propsychotic_significant_targets_all.tsv")
pp_int_sig_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/propsychotic_significant_targets_int.tsv")
ap_sig_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/antipsychotic_significant_targets.tsv")

uniqueN(pp_all_sig_targs) # [1] 172
uniqueN(pp_int_sig_targs) # [1] 117
uniqueN(ap_sig_targs) # [1] 129

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - DrugBank mechanism data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mechanisms <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/targets/drugbank_seachange_targets_with_actions.tsv")
mechanisms <- unique(mechanisms[,.(rxcui,drug_name,gene,action)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONVERT - Mechanism data into upregular or downregulator and remove everything else
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Downregualtion
mechanisms[action=="inhibitor" | 
		action=="antagonist" | 
		action=="blocker" | 
		action=="negative modulator" | 
		action=="inactivator" | 
		action=="suppressor" | 
		action=="weak inhibitor" | 
		action=="inhibitory allosteric modulator", action := "downregulator"]

# Upregulation
mechanisms[action=="agonist" | 
		action=="activator" | 
		action=="inducer" | 
		action=="potentiator" | 
		action=="positive allosteric modulator" | 
		action=="positive modulator" | 
		action=="stimulator", action := "upregulator"]

# Limit to only up/downregulation
mechanisms <- unique(mechanisms[action == "downregulator" | action=="upregulator",])

# Add combined symbol and mechanism column
mechanisms[,gene_mech := paste(gene,action,sep="_")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LIMIT - Mechanism data to significant targets only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mechanisms_pp_all <- unique(mechanisms[gene %in% pp_all_sig_targs$gene])
mechanisms_pp_int <- unique(mechanisms[gene %in% pp_int_sig_targs$gene])
mechanisms_ap <- unique(mechanisms[gene %in% ap_sig_targs$gene])

uniqueN(mechanisms_pp_all$gene_mech) # [1] 220
uniqueN(mechanisms_pp_int$gene_mech) # [1] 158
uniqueN(mechanisms_ap$gene_mech) # [1] 134

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load APs and PPs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_all <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_propsychotics_drugs_only.tsv")
pp_int <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/propsychotics_final.tsv")
ap <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/antipsychotics/ap_drugs.tsv")

pp_all[,rxcui := as.character(rxcui)]
pp_int[,rxcui := as.character(rxcui)]
ap[,rxcui := as.character(rxcui)]

uniqueN(pp_all$rxcui) # [1] 276
uniqueN(pp_int$rxcui) # [1] 147
uniqueN(ap$rxcui) # [1] 64

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load PP random sampling (all)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load
vgb <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_single_drug_filter_with_targets_drugs_only.tsv")
vgb[,rxcui := as.character(rxcui)]

# Limit to only drugs with mechanism data
vgb <- vgb[rxcui %in% mechanisms$rxcui,]
vgb[!rxcui %in% mechanisms$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name

# Remove all antipsychotic drugs
vgb <- vgb[!rxcui %in% ap$rxcui,]

# Save
pp_all_random_sampling <- unique(vgb$rxcui)
length(pp_all_random_sampling) # [1] 1635

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load PP random sampling (intersect)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load
sdr <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/sider_original_with_targets.tsv")
sdr[,rxcui := as.character(rxcui)]

# Limit to only drugs with mechanism data
sdr <- sdr[rxcui %in% mechanisms$rxcui,]
sdr[!rxcui %in% mechanisms$rxcui,] # Empty data.table (0 rows and 4 cols): rxcui,meddra,meddra_level,meddra_types

# Remove all antipsychotic drugs
sdr <- sdr[!rxcui %in% ap$rxcui,]
sdr_drugs <- unique(sdr$rxcui)
length(sdr_drugs) # [1] 1022

# Save
pp_int_random_sampling <- unique(intersect(pp_all_random_sampling,sdr_drugs))
pp_int_random_sampling <- unique(intersect(pp_int_random_sampling,mechanisms$rxcui))
length(pp_int_random_sampling) # [1] 1009

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load AP random sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load
atc <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/antipsychotics/atc_rxnorm_mapped.tsv")

# Limit to only drugs with mechanism data
atc_w_mechs <- atc[rxcui %in% mechanisms$rxcui,]

# Remove all propsychotic drugs
atc_w_mechs <- atc_w_mechs[!rxcui %in% pp_all$rxcui & !rxcui %in% pp_int$rxcui,]

# Save
ap_random_sampling <- unique(atc_w_mechs$rxcui)
length(ap_random_sampling) # [1] 1190

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Limit APs and PPs to significant targets with mechanism data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_all <- pp_all[rxcui %in% mechanisms_pp_all$rxcui,]
pp_int <- pp_int[rxcui %in% mechanisms_pp_int$rxcui,]
ap <- ap[rxcui %in% mechanisms_ap$rxcui,]

uniqueN(pp_all$rxcui) # [1] 167
uniqueN(pp_int$rxcui) # [1] 107
uniqueN(ap$rxcui) # [1] 44

pp_all[! rxcui %in% mechanisms_pp_all$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name
pp_int[! rxcui %in% mechanisms_pp_int$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name
ap[! rxcui %in% mechanisms_ap$rxcui,] # Empty data.table (0 rows and 2 cols): rxcui,drug_name

num_pp_all_observed <- uniqueN(pp_all$rxcui)
num_pp_int_observed <- uniqueN(pp_int$rxcui)
num_ap_observed <- uniqueN(ap$rxcui)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LIMIT - Mechanism data to drug targets only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mechanisms_pp_all <- unique(mechanisms_pp_all[rxcui %in% pp_all$rxcui])
mechanisms_pp_int <- unique(mechanisms_pp_int[rxcui %in% pp_int$rxcui])
mechanisms_ap <- unique(mechanisms_ap[rxcui %in% ap$rxcui])

uniqueN(mechanisms_pp_all$gene_mech) # [1] 164
uniqueN(mechanisms_pp_int$gene_mech) # [1] 116
uniqueN(mechanisms_ap$gene_mech) # [1] 51

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set number of permutations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perms <- 100000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE - PP all observed proportions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up the matrix 
mymtx_pp_all <- matrix(0, nrow=length(unique(mechanisms_pp_all$gene_mech)), ncol=(perms+3))
rownames(mymtx_pp_all) <- unique(mechanisms_pp_all$gene_mech)

# Add observed proportions to the first column
cur <- mechanisms_pp_all[rxcui %in% unique(pp_all$rxcui)][,list(frac=.N/num_pp_all_observed),gene_mech]
mymtx_pp_all[cur$gene_mech,1] <- cur$frac

# Add random sample proportions
for (i in 1:perms) {
  pp_all_random_set <- sample(pp_all_random_sampling, num_pp_all_observed, replace=FALSE, prob=NULL)
  cur <- mechanisms_pp_all[rxcui %in% pp_all_random_set][,list(frac=.N/num_pp_all_observed),gene_mech]
  mymtx_pp_all[cur$gene_mech,(i+1)] <- cur$frac
}

# Add mean random proportion for each gene
mymtx_pp_all[,(perms+2)] <- apply(mymtx_pp_all[,2:(perms+1)], 1, mean)

# Calculate p-value for each gene
mymtx_pp_all[,(perms+3)] <- apply(mymtx_pp_all, 1, function (x) length(which(x[2:(perms+1)] >= x[1]))/(perms))

# Check
mymtx_pp_all[1:5,c(1,(perms+2),(perms+3))]

# Save
save(mymtx_pp_all, file="/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/pp_all_significant_mechanisms.Rdata")

# Check that it loads
# load("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/pp_all_significant_mechanisms.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE - PP int observed proportions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up the matrix
mymtx_pp_int <- matrix(0, nrow=length(unique(mechanisms_pp_int$gene_mech)), ncol=(perms+3))
rownames(mymtx_pp_int) <- unique(mechanisms_pp_int$gene_mech)

# Add observed proportions to the first column
cur <- mechanisms_pp_int[rxcui %in% unique(pp_int$rxcui)][,list(frac=.N/num_pp_int_observed),gene_mech]
mymtx_pp_int[cur$gene_mech,1] <- cur$frac

# Add random sample proportions
for (i in 1:perms) {
  pp_int_random_set <- sample(pp_int_random_sampling, num_pp_int_observed, replace=FALSE, prob=NULL)
  cur <- mechanisms_pp_int[rxcui %in% pp_int_random_set][,list(frac=.N/num_pp_int_observed),gene_mech]
  mymtx_pp_int[cur$gene_mech,(i+1)] <- cur$frac
}

# Add mean random proportion for each gene
mymtx_pp_int[,(perms+2)] <- apply(mymtx_pp_int[,2:(perms+1)], 1, mean)

# Calculate p-value for each gene
mymtx_pp_int[,(perms+3)] <- apply(mymtx_pp_int, 1, function (x) length(which(x[2:(perms+1)] >= x[1]))/(perms))

# Check
mymtx_pp_int[1:5,c(1,(perms+2),(perms+3))]

# Save
save(mymtx_pp_int, file="/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/pp_int_significant_mechanisms.Rdata")

# Check that it loads
# load("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/pp_int_significant_mechanisms.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE - AP observed proportions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up the matrix
mymtx_ap <- matrix(0, nrow=length(unique(mechanisms_ap$gene_mech)), ncol=(perms+3))
rownames(mymtx_ap) <- unique(mechanisms_ap$gene_mech)

# Add observed proportions to the first column
cur <- mechanisms_ap[rxcui %in% unique(ap$rxcui)][,list(frac=.N/num_ap_observed),gene_mech]
mymtx_ap[cur$gene_mech,1] <- cur$frac

# Add random sample proportions
for (i in 1:perms) {
  ap_random_set <- sample(ap_random_sampling, num_ap_observed, replace=FALSE, prob=NULL)
  cur <- mechanisms_ap[rxcui %in% ap_random_set][,list(frac=.N/num_ap_observed),gene_mech]
  mymtx_ap[cur$gene_mech,(i+1)] <- cur$frac
}

# Add mean random proportion for each gene
mymtx_ap[,(perms+2)] <- apply(mymtx_ap[,2:(perms+1)], 1, mean)

# Calculate p-value for each gene
mymtx_ap[,(perms+3)] <- apply(mymtx_ap, 1, function (x) length(which(x[2:(perms+1)] >= x[1]))/(perms))

# Check
mymtx_ap[1:5,c(1,(perms+2),(perms+3))]

# Save
save(mymtx_ap, file="/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/ap_significant_mechanisms.Rdata")

# Check that it loads
# load("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/ap_significant_mechanisms.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract relevant data and convert to data tables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pp_all_sig_mechs <- data.table(gene_mech=rownames(mymtx_pp_all),
								observed_prop=mymtx_pp_all[,1],
								mean_rand_prop=mymtx_pp_all[,(perms+2)],
								p_value=mymtx_pp_all[,(perms+3)])

pp_int_sig_mechs <- data.table(gene_mech=rownames(mymtx_pp_int),
								observed_prop=mymtx_pp_int[,1],
								mean_rand_prop=mymtx_pp_int[,(perms+2)],
								p_value=mymtx_pp_int[,(perms+3)])

ap_sig_mechs <- data.table(gene_mech=rownames(mymtx_ap),
								observed_prop=mymtx_ap[,1],
								mean_rand_prop=mymtx_ap[,(perms+2)],
								p_value=mymtx_ap[,(perms+3)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bonferroni correction
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Establish corrected p-value
bonferroni_corrected_p <- (0.05/uniqueN(mechanisms$gene_mech))
bonferroni_corrected_p # [1] 2.061856e-05

# Save
pp_all_sig_mechs_bonferroni_correction <- pp_all_sig_mechs[p_value < bonferroni_corrected_p,]
pp_int_sig_mechs_bonferroni_correction <- pp_int_sig_mechs[p_value < bonferroni_corrected_p,]
ap_sig_mechs_bonferroni_correction <- ap_sig_mechs[p_value < bonferroni_corrected_p,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Isolate targets with p-value < 0.05
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_all_sig_mechs_p_05 <- pp_all_sig_mechs[p_value < 0.05,]
pp_int_sig_mechs_p_05 <- pp_int_sig_mechs[p_value < 0.05,]
ap_sig_mechs_p_05 <- ap_sig_mechs[p_value < 0.05,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check counts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check counts at p-value less than 0.05
uniqueN(pp_all_sig_mechs_p_05$gene) # [1] 115
uniqueN(pp_int_sig_mechs_p_05$gene) # [1] 89
uniqueN(ap_sig_mechs_p_05$gene) # [1] 51

# Check counts at corrected p-value
uniqueN(pp_all_sig_mechs[p_value < bonferroni_corrected_p,]) # [1] 74
uniqueN(pp_int_sig_mechs[p_value < bonferroni_corrected_p,]) # [1] 63
uniqueN(ap_sig_mechs[p_value < bonferroni_corrected_p,]) # [1] 29

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(pp_all_sig_mechs, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/propsychotic_significant_mechanisms_all.tsv")
fwrite(pp_int_sig_mechs, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/propsychotic_significant_mechanisms_int.tsv")
fwrite(ap_sig_mechs, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_mechanisms/antipsychotic_significant_mechanisms.tsv")
