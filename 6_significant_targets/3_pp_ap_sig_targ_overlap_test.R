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
perms <- 100000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load significant target results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/ap_significant_targets.Rdata")
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_all_significant_targets.Rdata")
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_int_significant_targets.Rdata")
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_no_cns_significant_targets.Rdata")

pp_all_sig_targs <- data.table(gene=rownames(mymtx_pp_all),
								observed_prop=mymtx_pp_all[,1],
								mean_rand_prop=mymtx_pp_all[,(perms+2)],
								p_value=mymtx_pp_all[,(perms+3)])

pp_int_sig_targs <- data.table(gene=rownames(mymtx_pp_int),
								observed_prop=mymtx_pp_int[,1],
								mean_rand_prop=mymtx_pp_int[,(perms+2)],
								p_value=mymtx_pp_int[,(perms+3)])

ap_sig_targs <- data.table(gene=rownames(mymtx_ap),
								observed_prop=mymtx_ap[,1],
								mean_rand_prop=mymtx_ap[,(perms+2)],
								p_value=mymtx_ap[,(perms+3)])

pp_no_cns_sig_targs <- data.table(gene=rownames(mymtx_pp_no_cns),
								observed_prop=mymtx_pp_no_cns[,1],
								mean_rand_prop=mymtx_pp_no_cns[,(perms+2)],
								p_value=mymtx_pp_no_cns[,(perms+3)])

uniqueN(pp_all_sig_targs[p_value < 0.05,]) # [1] 172
uniqueN(pp_int_sig_targs[p_value < 0.05,]) # [1] 117
uniqueN(ap_sig_targs[p_value < 0.05,]) # [1] 129
uniqueN(pp_no_cns_sig_targs[p_value < 0.05,]) # [1] 118

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create analysis table
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Isolate genes hit by either drug set
genes_all_pp <- unique(union(pp_all_sig_targs$gene,ap_sig_targs$gene))
genes_no_cns_pp <- unique(union(pp_no_cns_sig_targs$gene,ap_sig_targs$gene)) # Limit this to only targets of all non-CNS drugs

# Create table with main column for genes
myres_all_pp <- data.table(gene=genes_all_pp,pp_all_sig_targ=0,ap_sig_targ=0)
myres_no_cns_pp <- data.table(gene=genes_no_cns_pp,pp_no_cns_sig_targ=0,ap_sig_targ=0)

# Assign 1's to any gene's in a particular signifiant target set
myres_all_pp[gene %in% pp_all_sig_targs[p_value < 0.05,]$gene, pp_all_sig_targ := 1]
myres_all_pp[gene %in% ap_sig_targs[p_value < 0.05,]$gene, ap_sig_targ := 1]
myres_no_cns_pp[gene %in% pp_no_cns_sig_targs[p_value < 0.05,]$gene, pp_no_cns_sig_targ := 1]
myres_no_cns_pp[gene %in% ap_sig_targs[p_value < 0.05,]$gene, ap_sig_targ := 1]

uniqueN(myres_all_pp[pp_all_sig_targ==1,]$gene) # [1] 172
uniqueN(myres_all_pp[ap_sig_targ==1,]$gene) # [1] 129
uniqueN(myres_no_cns_pp[pp_no_cns_sig_targ==1,]$gene) # [1] 118
uniqueN(myres_no_cns_pp[ap_sig_targ==1,]$gene) # [1] 129

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Fisher's tests
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fisher's test for all PP significant targets
fishers_test_all_pp <- fisher.test(table(myres_all_pp$pp_all_sig_targ, myres_all_pp$ap_sig_targ))
fishers_p_all_pp <- fishers_test_all_pp$p.value
fishers_or_all_pp <- fishers_test_all_pp$estimate
add_all_pp <- data.table(target_sets="all_pp_vs_ap_sig_targs", odds_ratio=fishers_or_all_pp, p_value=fishers_p_all_pp)

# Fisher's test for no CNS PP significant targets
fishers_test_no_cns_pp <- fisher.test(table(myres_no_cns_pp$pp_no_cns_sig_targ, myres_no_cns_pp$ap_sig_targ))
fishers_p_no_cns_pp <- fishers_test_no_cns_pp$p.value
fishers_or_no_cns_pp <- fishers_test_no_cns_pp$estimate
add_no_cns_pp <- data.table(target_sets="no_cns_pp_vs_ap_sig_targs", odds_ratio=fishers_or_no_cns_pp, p_value=fishers_p_no_cns_pp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final <- rbind(add_all_pp,add_no_cns_pp)
fwrite(final, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_ap_overlap_significance.tsv")






