# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Isolate targets of non-CNS drugs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- unique(fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/targets/drugbank_seachange_targets_final.tsv"))
tar[,rxcui := as.character(rxcui)]

no_cns_drugs <- unique(fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/antipsychotics/atc_rxnorm_mapped.tsv")[,.(rxcui,drug_name,level1_description)])
no_cns_drugs <- no_cns_drugs[level1_description != "NERVOUS SYSTEM",]
no_cns_drugs <- unique(no_cns_drugs[,.(rxcui,drug_name)])
no_cns_drugs[,rxcui := as.character(rxcui)]

no_cns_targs <- unique(tar[rxcui %in% no_cns_drugs$rxcui]$gene)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load significant target results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp_all_sig_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/propsychotic_significant_targets_all.tsv")
pp_int_sig_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/propsychotic_significant_targets_int.tsv")
ap_sig_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/antipsychotic_significant_targets.tsv")
pp_no_cns_sig_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/propsychotic_significant_targets_no_cns.tsv")
pp_all_sig_targs_drugbank_only <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/propsychotic_significant_targets_all_drugbank_only.tsv")
pp_int_sig_targs_drugbank_only <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/propsychotic_significant_targets_int_drugbank_only.tsv")
ap_sig_targs_drugbank_only <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/antipsychotic_significant_targets_drugbank_only.tsv")
pp_no_cns_sig_targs_drugbank_only <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/propsychotic_significant_targets_no_cns_drugbank_only.tsv")

# Limit the non-CNS gene lists to target of non-CNS drugs
pp_no_cns_sig_targs <- pp_no_cns_sig_targs[gene %in% no_cns_targs,]
pp_no_cns_sig_targs_drugbank_only <- pp_no_cns_sig_targs_drugbank_only[gene %in% no_cns_targs,]

# Check counts - all
uniqueN(pp_all_sig_targs) # [1] 3402
uniqueN(pp_int_sig_targs) # [1] 3402
uniqueN(ap_sig_targs) # [1] 3402
uniqueN(pp_no_cns_sig_targs) # [1] 2495
uniqueN(pp_all_sig_targs_drugbank_only) # [1] 2760
uniqueN(pp_int_sig_targs_drugbank_only) # [1] 2760
uniqueN(ap_sig_targs_drugbank_only) # [1] 2760
uniqueN(pp_no_cns_sig_targs_drugbank_only) # [1] 1861

# Check counts - significant
uniqueN(pp_all_sig_targs[p_value < 0.05,]) # [1] 169
uniqueN(pp_int_sig_targs[p_value < 0.05,]) # [1] 118
uniqueN(ap_sig_targs[p_value < 0.05,]) # [1] 129
uniqueN(pp_no_cns_sig_targs[p_value < 0.05,]) # [1] 118
uniqueN(pp_all_sig_targs_drugbank_only[p_value < 0.05,]) # [1] 86
uniqueN(pp_int_sig_targs_drugbank_only[p_value < 0.05,]) # [1] 78
uniqueN(ap_sig_targs_drugbank_only[p_value < 0.05,]) # [1] 38
uniqueN(pp_no_cns_sig_targs_drugbank_only[p_value < 0.05,]) # [1] 44

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create analysis table - DrugBank and SeaChange
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Isolate genes hit by either drug set
genes_all_pp <- unique(intersect(pp_all_sig_targs$gene,ap_sig_targs$gene))
genes_no_cns_pp <- unique(intersect(pp_no_cns_sig_targs$gene,ap_sig_targs$gene))

length(genes_all_pp) # [1] 3402
length(genes_no_cns_pp) # [1] 2495

# Create table with main column for genes
myres_all_pp <- data.table(gene=genes_all_pp,pp_all_sig_targ=0,ap_sig_targ=0)
myres_no_cns_pp <- data.table(gene=genes_no_cns_pp,pp_no_cns_sig_targ=0,ap_sig_targ=0)

# Assign 1's to any gene's in a particular signifiant target set
myres_all_pp[gene %in% pp_all_sig_targs[p_value < 0.05,]$gene, pp_all_sig_targ := 1]
myres_all_pp[gene %in% ap_sig_targs[p_value < 0.05,]$gene, ap_sig_targ := 1]
myres_no_cns_pp[gene %in% pp_no_cns_sig_targs[p_value < 0.05,]$gene, pp_no_cns_sig_targ := 1]
myres_no_cns_pp[gene %in% ap_sig_targs[p_value < 0.05,]$gene, ap_sig_targ := 1]

uniqueN(myres_all_pp[pp_all_sig_targ==1,]$gene) # [1] 169
uniqueN(myres_all_pp[ap_sig_targ==1,]$gene) # [1] 129
uniqueN(myres_no_cns_pp[pp_no_cns_sig_targ==1,]$gene) # [1] 118
uniqueN(myres_no_cns_pp[ap_sig_targ==1,]$gene) # [1] 127

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create analysis table - DrugBank only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Isolate genes hit by either drug set
genes_all_pp_drugbank_only <- unique(intersect(pp_all_sig_targs_drugbank_only$gene,ap_sig_targs_drugbank_only$gene))
genes_no_cns_pp_drugbank_only <- unique(intersect(pp_no_cns_sig_targs_drugbank_only$gene,ap_sig_targs_drugbank_only$gene))

length(genes_all_pp_drugbank_only) # [1] 2760
length(genes_no_cns_pp_drugbank_only) # [1] 1861

# Create table with main column for genes
myres_all_pp_drugbank_only <- data.table(gene=genes_all_pp_drugbank_only,pp_all_sig_targ_drugbank_only=0,ap_sig_targ_drugbank_only=0)
myres_no_cns_pp_drugbank_only <- data.table(gene=genes_no_cns_pp_drugbank_only,pp_no_cns_sig_targ_drugbank_only=0,ap_sig_targ_drugbank_only=0)

# Assign 1's to any gene's in a particular signifiant target set
myres_all_pp_drugbank_only[gene %in% pp_all_sig_targs_drugbank_only[p_value < 0.05,]$gene, pp_all_sig_targ_drugbank_only := 1]
myres_all_pp_drugbank_only[gene %in% ap_sig_targs_drugbank_only[p_value < 0.05,]$gene, ap_sig_targ_drugbank_only := 1]
myres_no_cns_pp_drugbank_only[gene %in% pp_no_cns_sig_targs_drugbank_only[p_value < 0.05,]$gene, pp_no_cns_sig_targ_drugbank_only := 1]
myres_no_cns_pp_drugbank_only[gene %in% ap_sig_targs_drugbank_only[p_value < 0.05,]$gene, ap_sig_targ_drugbank_only := 1]

uniqueN(myres_all_pp_drugbank_only[pp_all_sig_targ_drugbank_only==1,]$gene) # [1] 86
uniqueN(myres_all_pp_drugbank_only[ap_sig_targ_drugbank_only==1,]$gene) # [1] 38
uniqueN(myres_no_cns_pp_drugbank_only[pp_no_cns_sig_targ_drugbank_only==1,]$gene) # [1] 44
uniqueN(myres_no_cns_pp_drugbank_only[ap_sig_targ_drugbank_only==1,]$gene) # [1] 36

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Fisher's tests - DrugBank and SeaChange
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
# Run Fisher's tests - DrugBank only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fisher's test for all PP significant targets
fishers_test_all_pp_drugbank_only <- fisher.test(table(myres_all_pp_drugbank_only$pp_all_sig_targ_drugbank_only, myres_all_pp_drugbank_only$ap_sig_targ_drugbank_only))
fishers_p_all_pp_drugbank_only <- fishers_test_all_pp_drugbank_only$p.value
fishers_or_all_pp_drugbank_only <- fishers_test_all_pp_drugbank_only$estimate
add_all_pp_drugbank_only <- data.table(target_sets="all_pp_vs_ap_sig_targs_drugbank_only", odds_ratio=fishers_or_all_pp_drugbank_only, p_value=fishers_p_all_pp_drugbank_only)

# Fisher's test for no CNS PP significant targets
fishers_test_no_cns_pp_drugbank_only <- fisher.test(table(myres_no_cns_pp_drugbank_only$pp_no_cns_sig_targ_drugbank_only, myres_no_cns_pp_drugbank_only$ap_sig_targ_drugbank_only))
fishers_p_no_cns_pp_drugbank_only <- fishers_test_no_cns_pp_drugbank_only$p.value
fishers_or_no_cns_pp_drugbank_only <- fishers_test_no_cns_pp_drugbank_only$estimate
add_no_cns_pp_drugbank_only <- data.table(target_sets="no_cns_pp_vs_ap_sig_targs_drugbank_only", odds_ratio=fishers_or_no_cns_pp_drugbank_only, p_value=fishers_p_no_cns_pp_drugbank_only)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final <- rbind(add_all_pp,add_no_cns_pp)
final <- rbind(final,add_all_pp_drugbank_only)
final <- rbind(final,add_no_cns_pp_drugbank_only)
fwrite(final, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/significant_targets/pp_ap_overlap_significance_drugbank_only.tsv")






