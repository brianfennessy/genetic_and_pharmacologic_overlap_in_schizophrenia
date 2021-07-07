# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Store user home directory and navigate to it
user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
setwd(user_home_directory)

# Get the official home directory of the user
home_dir <- getwd()

# Navigate to the raw data
setwd(home_dir)

source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
library("weights")
options(width=200)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load significant targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_qc_sig_targs <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_all.tsv")
pp_qc_plus_sig_targs <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_int.tsv")
ap_sig_targs <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/antipsychotic_significant_targets.tsv")
pp_no_cns_sig_targs <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_no_cns.tsv")
pp_qc_sig_targs_drugbank_only <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_all_drugbank_only.tsv")
pp_qc_plus_sig_targs_drugbank_only <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_int_drugbank_only.tsv")
ap_sig_targs_drugbank_only <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/antipsychotic_significant_targets_drugbank_only.tsv")
pp_no_cns_sig_targs_drugbank_only <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_no_cns_drugbank_only.tsv")
all_drug_targs <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")
drugbank_targs <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")
drugbank_targs <- unique(drugbank_targs[database == "drugbank",])

uniqueN(pp_qc_sig_targs$gene) # [1] 3402
uniqueN(pp_qc_plus_sig_targs$gene) # [1] 3402
uniqueN(ap_sig_targs$gene) # [1] 3402
uniqueN(pp_no_cns_sig_targs$gene) # [1] 3402
uniqueN(pp_qc_sig_targs_drugbank_only$gene) # [1] 2760
uniqueN(pp_qc_plus_sig_targs_drugbank_only$gene) # [1] 2760
uniqueN(ap_sig_targs_drugbank_only$gene) # [1] 2760
uniqueN(pp_no_cns_sig_targs_drugbank_only$gene) # [1] 2760
uniqueN(all_drug_targs$gene) # [1] 3402
uniqueN(drugbank_targs$gene) # [1] 2760

pp_qc_targs <- unique(pp_qc_sig_targs[p_value < 0.05,]$gene)
pp_qc_plus_targs <- unique(pp_qc_plus_sig_targs[p_value < 0.05,]$gene)
ap_targs <- unique(ap_sig_targs[p_value < 0.05,]$gene)
pp_no_cns_targs <- unique(pp_no_cns_sig_targs[p_value < 0.05,]$gene)
pp_qc_targs_drugbank_only <- unique(pp_qc_sig_targs_drugbank_only[p_value < 0.05,]$gene)
pp_qc_plus_targs_drugbank_only <- unique(pp_qc_plus_sig_targs_drugbank_only[p_value < 0.05,]$gene)
ap_targs_drugbank_only <- unique(ap_sig_targs_drugbank_only[p_value < 0.05,]$gene)
pp_no_cns_targs_drugbank_only <- unique(pp_no_cns_sig_targs_drugbank_only[p_value < 0.05,]$gene)
all_drug_targs <- unique(all_drug_targs$gene)
drugbank_targs <- unique(drugbank_targs$gene)

length(pp_qc_targs) # [1] 169
length(pp_qc_plus_targs) # [1] 118
length(ap_targs) # [1] 129
length(pp_no_cns_targs) # [1] 118
length(pp_qc_targs_drugbank_only) # [1] 86
length(pp_qc_plus_targs_drugbank_only) # [1] 78
length(ap_targs_drugbank_only) # [1] 38
length(pp_no_cns_targs_drugbank_only) # [1] 44
length(all_drug_targs) # [1] 3402
length(drugbank_targs) # [1] 2760

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load SCHEMA, PrediXcan, and PGC nearest gene data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("/sc/arion/projects/psychgen/methods/rx/data/molPsych2019_psychosisGenData.Rdata") #genlist

# SCHEMA
sma <- genlist$sma
sma <- unique(sma[!is.na(gene_name)])

# MAGMA
mag <- genlist$mag
mag <- unique(mag[!is.na(GENE)])

# Drug Targets from DrugBank and SeaChange
drug_targets_for_nearest_gene <- data.table(symbol=all_drug_targs)

# PGC3 nearest gene
# pgc_nearest_gene <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/sz_targets/nearest_gene/nearest_gene_common_variants.tsv")
# setnames(pgc_nearest_gene,"gene","symbol")
# pgc_near <- copy(pgc_nearest_gene)
# pgc_near[,neglog10.pvalue := log10(p_value) * (-1)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initiate results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myres <- c()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initiate target sets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
targ_sets <- list(all_drug_targs, drugbank_targs)
names(targ_sets) <- c("all_drug_targs","drugbank_targs")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Testing for enrichment of all drug targets in SCHEMA
#
# Rows: Yes/no, is it a significant SCHEMA gene? (SIG == 1 / SIG == 0)
# Columns: Yes/no, is it a drug target in DrugBank or SeaChange? (TARG == 1 / TARG == 0)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# Calculate bonferroni threshold for SCHEMA
schema_bonferroni_corrected <- (0.05/uniqueN(sma$gene_name))

for (threshold in c(0.05,schema_bonferroni_corrected)) {

	for (set in names(targ_sets)) {

		# Save current set
		cur_set <- targ_sets[[set]]

		# Yes/no, is it a significant SCHEMA gene? (SIG == 1 / SIG == 0)
		sma[, SIG := 0]
		sma[p_studywide < threshold, SIG := 1]

		# Yes/no, is it a drug target in DrugBank or SeaChange? (TARG == 1 / TARG == 0)
		sma[, TARG := 0]
		sma[gene_name %in% cur_set, TARG := 1]

		# Wilcox
		wilcox_p <- wilcox.test(sma[TARG==1]$neglog10.p_studywide, sma[TARG==0]$neglog10.p_studywide, alternative="greater")$p.value

		# Weighted correlation test
		weighted_correlation_test <- wtd.cor(sma$TARG, sma$SIG)
		weighted_correlation <- weighted_correlation_test[1]
		weighted_correlation_p <- weighted_correlation_test[4]

		# Fisher's test
		fishers_test <- fisher.test(table(sma$TARG, sma$SIG))
		fishers_p <- fishers_test$p.value
		fishers_or <- fishers_test$estimate

		# Get counts for sanity check
		or_a = nrow(sma[TARG==1 & SIG==1])
		or_b = nrow(sma[TARG==1 & SIG==0])
		or_c = nrow(sma[TARG==0 & SIG==1])
		or_d = nrow(sma[TARG==0 & SIG==0])
		or = (or_a/or_b) / (or_c/or_d)

		# Add to results table
		add <- data.table(targlist=set, genset="schema", gene_sig_threshold=threshold, 
						wilcox_p=wilcox_p,
						weighted_correlation=weighted_correlation, weighted_correlation_p=weighted_correlation_p,
						fishers_or=fishers_or, fishers_p=fishers_p,
						targ1=or_a+or_b, sig1=or_a+or_c,
						targ1scz1=or_a, targ1scz0=or_b,
						targ0scz1=or_c, targ0scz0=or_d)
		myres <- rbind(myres, add)

	}

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Testing for enrichment of all drug targets in MAGMA
#
# Rows: Yes/no, is it a significant MAGMA gene? (SIG == 1 / SIG == 0)
# Columns: Yes/no, is it a drug target in DrugBank or SeaChange? (TARG == 1 / TARG == 0)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# Calculate bonferroni threshold for MAGMA
magma_bonferroni_corrected <- (0.05/uniqueN(mag$GENE))

for (threshold in c(0.05,magma_bonferroni_corrected)) {

	for (set in names(targ_sets)) {

		# Save current set
		cur_set <- targ_sets[[set]]

		# Yes/no, is it a significant MAGMA gene? (SIG == 1 / SIG == 0)
		mag[, SIG := 0]
		mag[P < threshold, SIG := 1]

		# Yes/no, is it a drug target in DrugBank or SeaChange? (TARG == 1 / TARG == 0)
		mag[, TARG := 0]
		mag[GENE %in% cur_set, TARG := 1]

		# Wilcox
		wilcox_p <- wilcox.test(mag[TARG==1]$neglog10.P, mag[TARG==0]$neglog10.P, alternative="greater")$p.value

		# Weighted correlation test
		weighted_correlation_test <- wtd.cor(mag$TARG, mag$SIG)
		weighted_correlation <- weighted_correlation_test[1]
		weighted_correlation_p <- weighted_correlation_test[4]

		# Fisher's at p < 0.05
		fishers_test <- fisher.test(table(mag$TARG, mag$SIG))
		fishers_p <- fishers_test$p.value
		fishers_or <- fishers_test$estimate

		# Get counts for sanity check
		or_a = nrow(mag[TARG==1 & SIG==1])
		or_b = nrow(mag[TARG==1 & SIG==0])
		or_c = nrow(mag[TARG==0 & SIG==1])
		or_d = nrow(mag[TARG==0 & SIG==0])
		or = (or_a/or_b) / (or_c/or_d)

		# Add to results table
		add <- data.table(targlist=set, genset="magma", gene_sig_threshold=threshold,
							wilcox_p=wilcox_p,
							weighted_correlation=weighted_correlation, weighted_correlation_p=weighted_correlation_p,
							fishers_or=fishers_or, fishers_p=fishers_p,
							targ1=or_a+or_b, sig1=or_a+or_c,
							targ1scz1=or_a, targ1scz0=or_b,
							targ0scz1=or_c, targ0scz0=or_d)
		myres <- rbind(myres, add)

	}

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Limit SCHEMA and MAGMA to symbols in DrugBank or SeaChange for further testing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sma <- sma[gene_name %in% all_drug_targs,]
mag <- mag[GENE %in% all_drug_targs,]

targ_sets <- list(pp_qc_targs,pp_qc_plus_targs,ap_targs,pp_no_cns_targs,
					pp_qc_targs_drugbank_only,pp_qc_plus_targs_drugbank_only,ap_targs_drugbank_only,pp_no_cns_targs_drugbank_only)

names(targ_sets) <- c("pp_qc_targs","pp_qc_plus_targs","ap_targs","pp_no_cns_targs",
						"pp_qc_targs_drugbank_only","pp_qc_plus_targs_drugbank_only","ap_targs_drugbank_only","pp_no_cns_targs_drugbank_only")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Testing for enrichment of QC PP / QC+ PP / AP targets in SCHEMA (limited to DrugBank/SeaChange)
#
# Rows: Yes/no, is it a significant SCHEMA gene? (SIG == 1 / SIG == 0)
# Columns: Yes/no, is it a significant drug target? (TARG == 1 / TARG == 0)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

for (threshold in c(0.05,schema_bonferroni_corrected)) {

	for (set in names(targ_sets)) {

		# Save current set
		cur_set <- targ_sets[[set]]

		# Yes/no, is it a significant SCHEMA gene? (SIG == 1 / SIG == 0)
		sma[, SIG := 0]
		sma[p_studywide < threshold, SIG := 1]

		# Yes/no, is it a significant drug target? (TARG == 1 / TARG == 0)
		sma[, TARG := 0]
		sma[gene_name %in% cur_set, TARG := 1]

		# Wilcox
		wilcox_p <- wilcox.test(sma[TARG==1]$neglog10.p_studywide, sma[TARG==0]$neglog10.p_studywide, alternative="greater")$p.value

		# Weighted correlation test
		weighted_correlation_test <- wtd.cor(sma$TARG, sma$SIG)
		weighted_correlation <- weighted_correlation_test[1]
		weighted_correlation_p <- weighted_correlation_test[4]

		# Fisher's at p < 0.05
		fishers_test <- fisher.test(table(sma$TARG, sma$SIG))
		fishers_p <- fishers_test$p.value
		fishers_or <- fishers_test$estimate

		# Get counts for sanity check
		or_a = nrow(sma[TARG==1 & SIG==1])
		or_b = nrow(sma[TARG==1 & SIG==0])
		or_c = nrow(sma[TARG==0 & SIG==1])
		or_d = nrow(sma[TARG==0 & SIG==0])
		or = (or_a/or_b) / (or_c/or_d)

		# Add to results table
		add <- data.table(targlist=set, genset="schema", gene_sig_threshold=threshold,
							wilcox_p=wilcox_p,
							weighted_correlation=weighted_correlation, weighted_correlation_p=weighted_correlation_p,
							fishers_or=fishers_or, fishers_p=fishers_p,
							targ1=or_a+or_b, sig1=or_a+or_c,
							targ1scz1=or_a, targ1scz0=or_b,
							targ0scz1=or_c, targ0scz0=or_d)
		myres <- rbind(myres, add)

	}

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Testing for enrichment of QC PP / QC+ PP / AP targets in MAGMA (limited to DrugBank/SeaChange)
#
# Rows: Yes/no, is it a significant MAGMA gene? (SIG == 1 / SIG == 0)
# Columns: Yes/no, is it a drug target? (TARG == 1 / TARG == 0)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

for (threshold in c(0.05,magma_bonferroni_corrected)) {

	for (set in names(targ_sets)) {

		# Save current set
		cur_set <- targ_sets[[set]]

		# Yes/no, is it a significant MAGMA gene? (SIG == 1 / SIG == 0)
		mag[, SIG := 0]
		mag[P < threshold, SIG := 1]

		# Yes/no, is it a significant drug target? (TARG == 1 / TARG == 0)
		mag[, TARG := 0]
		mag[GENE %in% cur_set, TARG := 1]

		# Wilcox
		wilcox_p <- wilcox.test(mag[TARG==1]$neglog10.P, mag[TARG==0]$neglog10.P, alternative="greater")$p.value

		# Weighted correlation test
		weighted_correlation_test <- wtd.cor(mag$TARG, mag$SIG)
		weighted_correlation <- weighted_correlation_test[1]
		weighted_correlation_p <- weighted_correlation_test[4]

		# Fisher's at p < 0.05
		fishers_test <- fisher.test(table(mag$TARG, mag$SIG))
		fishers_p <- fishers_test$p.value
		fishers_or <- fishers_test$estimate

		# Get counts for sanity check
		or_a = nrow(mag[TARG==1 & SIG==1])
		or_b = nrow(mag[TARG==1 & SIG==0])
		or_c = nrow(mag[TARG==0 & SIG==1])
		or_d = nrow(mag[TARG==0 & SIG==0])
		or = (or_a/or_b) / (or_c/or_d)

		# Add to results table
		add <- data.table(targlist=set, genset="magma", gene_sig_threshold=threshold,
							wilcox_p=wilcox_p,
							weighted_correlation=weighted_correlation, weighted_correlation_p=weighted_correlation_p,
							fishers_or=fishers_or, fishers_p=fishers_p,
							targ1=or_a+or_b, sig1=or_a+or_c,
							targ1scz1=or_a, targ1scz0=or_b,
							targ0scz1=or_c, targ0scz0=or_d)
		myres <- rbind(myres, add)

	}

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Testing for enrichment of QC PP / QC+ PP / AP targets & PGC3 nearest gene in DrugBank/SeaChange
#
# Rows: Yes/no, is it a significant PGC3 nearest gene? (SIG == 1 / SIG == 0)
# Columns: Yes/no, is it a significant drug target? (TARG == 1 / TARG == 0)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# for (set in names(targ_sets)) {

# 	# Save current set
# 	cur_set <- targ_sets[[set]]

# 	# Yes/no, is it a significant PGC3 nearest gene? (SIG == 1 / SIG == 0)
# 	drug_targets_for_nearest_gene[,SIG:=0]
#     drug_targets_for_nearest_gene[symbol %in% pgc_near$symbol,SIG:=1]

# 	# Yes/no, is it a significant drug target? (TARG == 1 / TARG == 0)
# 	drug_targets_for_nearest_gene[, TARG := 0]
# 	drug_targets_for_nearest_gene[symbol %in% cur_set, TARG := 1]

# 	# Wilcox
# 	wilcox_p <- "Not Applicable"

# 	# Weighted correlation test
# 	weighted_correlation_test <- wtd.cor(drug_targets_for_nearest_gene$TARG, drug_targets_for_nearest_gene$SIG)
# 	weighted_correlation <- weighted_correlation_test[1]
# 	weighted_correlation_p <- weighted_correlation_test[4]

# 	# Fisher's at p < 0.05
# 	fishers_test <- fisher.test(table(drug_targets_for_nearest_gene$TARG, drug_targets_for_nearest_gene$SIG))
# 	fishers_p <- fishers_test$p.value
# 	fishers_or <- fishers_test$estimate

# 	# Get counts for sanity check
# 	or_a = nrow(drug_targets_for_nearest_gene[TARG==1 & SIG==1])
# 	or_b = nrow(drug_targets_for_nearest_gene[TARG==1 & SIG==0])
# 	or_c = nrow(drug_targets_for_nearest_gene[TARG==0 & SIG==1])
# 	or_d = nrow(drug_targets_for_nearest_gene[TARG==0 & SIG==0])
# 	or = (or_a/or_b) / (or_c/or_d)

# 	# Add to results table
# 	add <- data.table(targlist=set, genset="nearest_gene_to_pgc3_snp", gene_sig_threshold= "Not Applicable",
# 						wilcox_p=wilcox_p,
# 						weighted_correlation=weighted_correlation, weighted_correlation_p=weighted_correlation_p,
# 						fishers_or=fishers_or, fishers_p=fishers_p,
# 						targ1=or_a+or_b, sig1=or_a+or_c,
# 						targ1scz1=or_a, targ1scz0=or_b,
# 						targ0scz1=or_c, targ0scz0=or_d)
# 	myres <- rbind(myres, add)

# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Take a looksie!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setorder(myres,targlist)
myres[,.(targlist,genset,gene_sig_threshold,fishers_or,fishers_p,wilcox_p,weighted_correlation,weighted_correlation_p)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(myres, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/sz_targets/enrichment_results.tsv")







