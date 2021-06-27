# ==================================================================================================================================================================
#
# GET - Libraries
#
# ==================================================================================================================================================================

# Store user home directory and navigate to it
user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
setwd(user_home_directory)

# Get the official home directory of the user
home_dir <- getwd()

# Navigate to the raw data zip file and unzip
setwd(home_dir)

source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
library(xlsx)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(stringr)
options(width=160)

# ==================================================================================================================================================================
#
# SHEETS 1 & 2 - QC PP & QC+ PP drugs, side effects, and level 3 classifications
#
# ==================================================================================================================================================================

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_med_qc <- unique(fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/qc_pp_drugs_meddra_descriptions.tsv")[,.(rxcui,drug_name,adverse_event=pt_name)])
vgb_med_qc_plus <- unique(fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/qc_plus_pp_drugs_meddra_descriptions.tsv")[,.(rxcui,drug_name,adverse_event=pt_name)])

vgb_atc_qc <- unique(fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/qc_pp_drugs_atc_descriptions.tsv")[,.(rxcui,drug_name,drug_type=level3_description)])
vgb_atc_qc_plus <- unique(fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/qc_plus_pp_drugs_atc_descriptions.tsv")[,.(rxcui,drug_name,drug_type=level3_description)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MedDRA - Collapse all terms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uniqueN(vgb_med_qc$rxcui) # [1] 276
uniqueN(vgb_med_qc_plus$rxcui) # [1] 148

# Optional (counts if we don't want to collapse)
vgb_med_qc_counts <- vgb_med_qc[,.N,adverse_event]
vgb_med_qc_counts[,prop_qc := N/uniqueN(vgb_med_qc$rxcui)]
setorder(vgb_med_qc_counts,-prop_qc)

# Optional (counts if we don't want to collapse)
vgb_med_qc_plus_counts <- vgb_med_qc_plus[,.N,adverse_event]
vgb_med_qc_plus_counts[,prop_qc_plus := N/uniqueN(vgb_med_qc_plus$rxcui)]
setorder(vgb_med_qc_plus_counts,-prop_qc_plus)

# Collapse
vgb_med_qc <- vgb_med_qc[, lapply(.SD, paste0, collapse="|"), by=c("rxcui","drug_name")]
vgb_med_qc_plus <- vgb_med_qc_plus[, lapply(.SD, paste0, collapse="|"), by=c("rxcui","drug_name")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ATC - Collapse all terms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_atc_qc[,drug_type := str_to_sentence(drug_type)]
vgb_atc_qc_plus[,drug_type := str_to_sentence(drug_type)]

uniqueN(vgb_atc_qc$rxcui) # [1] 276
uniqueN(vgb_atc_qc_plus$rxcui) # [1] 148

# Optional (counts if we don't want to collapse)
vgb_atc_qc_counts <- vgb_atc_qc[,.N,drug_type]
vgb_atc_qc_counts[,prop_qc := N/uniqueN(vgb_atc_qc$rxcui)]
setorder(vgb_atc_qc_counts,-prop_qc)

# Optional (counts if we don't want to collapse)
vgb_atc_qc_plus_counts <- vgb_atc_qc_plus[,.N,drug_type]
vgb_atc_qc_plus_counts[,prop_qc_plus := N/uniqueN(vgb_atc_qc_plus$rxcui)]
setorder(vgb_atc_qc_plus_counts,-prop_qc_plus)

# Collapse
vgb_atc_qc <- vgb_atc_qc[, lapply(.SD, paste0, collapse="|"), by=c("rxcui","drug_name")]
vgb_atc_qc_plus <- vgb_atc_qc_plus[, lapply(.SD, paste0, collapse="|"), by=c("rxcui","drug_name")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge final tables - all information in one
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sheet1 <- merge(vgb_med_qc,vgb_atc_qc,by=c("rxcui","drug_name"))
sheet2 <- merge(vgb_med_qc_plus,vgb_atc_qc_plus,by=c("rxcui","drug_name"))

# ==================================================================================================================================================================
#
# SHEET 3 - AP drugs, side effects, and level 3 classifications
#
# ==================================================================================================================================================================
sheet3 <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/ap_drugs.tsv")

# ==================================================================================================================================================================
#
# SHEET 4 - Target proportions
#
# ==================================================================================================================================================================

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")
target_proportions <- unique(tar[,.(gene)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC Propsychotics targets (significant) - DrugBank and SeaChange
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qc_pp_targs_sig <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_all.tsv")[,.(gene,total_prop_db_and_sc=total_prop,observed_prop_qc_pp=observed_prop,p_value_qc_pp=p_value)]
uniqueN(qc_pp_targs_sig[p_value_qc_pp < 0.05,]$gene) # [1] 169

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC+ Propsychotics targets (significant) - DrugBank and SeaChange
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qc_plus_pp_targs_sig <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_int.tsv")[,.(gene,total_prop_db_and_sc=total_prop,observed_prop_qc_plus_pp=observed_prop,p_value_qc_plus_pp=p_value)]
uniqueN(qc_plus_pp_targs_sig[p_value_qc_plus_pp < 0.05,]$gene) # [1] 118

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Antipsychotics targets (significant) - DrugBank and SeaChange
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ap_targs_sig <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/antipsychotic_significant_targets.tsv")[,.(gene,total_prop_db_and_sc=total_prop,observed_prop_ap=observed_prop,p_value_ap=p_value)]
uniqueN(ap_targs_sig[p_value_ap < 0.05,]$gene) # [1] 129

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No CNS Propsychotics targets (significant) - DrugBank and SeaChange
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
no_cns_pp_targs_sig <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_no_cns.tsv")[,.(gene,total_prop_db_and_sc=total_prop,observed_prop_no_cns_pp=observed_prop,p_value_no_cns_pp=p_value)]
uniqueN(no_cns_pp_targs_sig[p_value_no_cns_pp < 0.05,]$gene) # [1] 118

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC Propsychotics targets (significant) - DrugBank only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qc_pp_targs_sig_drugbank_only <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_all_drugbank_only.tsv")[,.(gene,total_prop_db_only=total_prop,observed_prop_qc_pp_drugbank_only=observed_prop,p_value_qc_pp_drugbank_only=p_value)]
uniqueN(qc_pp_targs_sig_drugbank_only[p_value_qc_pp_drugbank_only < 0.05,]$gene) # [1] 86

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC+ Propsychotics targets (significant) - DrugBank only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qc_plus_pp_targs_sig_drugbank_only <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_int_drugbank_only.tsv")[,.(gene,total_prop_db_only=total_prop,observed_prop_qc_plus_pp_drugbank_only=observed_prop,p_value_qc_plus_pp_drugbank_only=p_value)]
uniqueN(qc_plus_pp_targs_sig_drugbank_only[p_value_qc_plus_pp_drugbank_only < 0.05,]$gene) # [1] 78

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Antipsychotics targets (significant) - DrugBank only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ap_targs_sig_drugbank_only <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/antipsychotic_significant_targets_drugbank_only.tsv")[,.(gene,total_prop_db_only=total_prop,observed_prop_ap_drugbank_only=observed_prop,p_value_ap_drugbank_only=p_value)]
uniqueN(ap_targs_sig_drugbank_only[p_value_ap_drugbank_only < 0.05,]$gene) # [1] 38

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No CNS Propsychotics targets (significant) - DrugBank only
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
no_cns_pp_targs_sig_drugbank_only <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/propsychotic_significant_targets_no_cns_drugbank_only.tsv")[,.(gene,total_prop_db_only=total_prop,observed_prop_no_cns_pp_drugbank_only=observed_prop,p_value_no_cns_pp_drugbank_only=p_value)]
uniqueN(no_cns_pp_targs_sig_drugbank_only[p_value_no_cns_pp_drugbank_only < 0.05,]$gene) # [1] 44

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_proportions <- merge(target_proportions,qc_pp_targs_sig,by="gene",all.x=TRUE)
target_proportions <- merge(target_proportions,qc_plus_pp_targs_sig,by=c("gene","total_prop_db_and_sc"),all.x=TRUE)
target_proportions <- merge(target_proportions,ap_targs_sig,by=c("gene","total_prop_db_and_sc"),all.x=TRUE)
target_proportions <- merge(target_proportions,no_cns_pp_targs_sig,by=c("gene","total_prop_db_and_sc"),all.x=TRUE)
target_proportions <- merge(target_proportions,qc_pp_targs_sig_drugbank_only,by="gene",all.x=TRUE)
target_proportions <- merge(target_proportions,qc_plus_pp_targs_sig_drugbank_only,by=c("gene","total_prop_db_only"),all.x=TRUE)
target_proportions <- merge(target_proportions,ap_targs_sig_drugbank_only,by=c("gene","total_prop_db_only"),all.x=TRUE)
target_proportions <- merge(target_proportions,no_cns_pp_targs_sig_drugbank_only,by=c("gene","total_prop_db_only"),all.x=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Relabel any NA's
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

for (i in c("total_prop_db_only", "total_prop_db_and_sc", 
				"observed_prop_qc_pp", "observed_prop_qc_plus_pp", "observed_prop_ap", "observed_prop_no_cns_pp",
				"observed_prop_qc_pp_drugbank_only", "observed_prop_qc_plus_pp_drugbank_only",
				"observed_prop_ap_drugbank_only", "observed_prop_no_cns_pp_drugbank_only")){target_proportions[is.na(get(i)), (i):=0]}

for (i in c("p_value_qc_pp", "p_value_qc_plus_pp", "p_value_ap", "p_value_no_cns_pp",
				"p_value_qc_pp_drugbank_only", "p_value_qc_plus_pp_drugbank_only",
				"p_value_ap_drugbank_only", "p_value_no_cns_pp_drugbank_only")){target_proportions[is.na(get(i)), (i):=1]}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sanity check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uniqueN(target_proportions[observed_prop_qc_pp != 0,]$gene) # [1] 1134
uniqueN(target_proportions[observed_prop_qc_pp != 0 & p_value_qc_pp < 0.05,]$gene) # [1] 169

uniqueN(target_proportions[observed_prop_qc_plus_pp != 0,]$gene) # [1] 963
uniqueN(target_proportions[observed_prop_qc_plus_pp != 0 & p_value_qc_plus_pp < 0.05,]$gene) # [1] 118

uniqueN(target_proportions[observed_prop_no_cns_pp != 0,]$gene) # [1] 916
uniqueN(target_proportions[observed_prop_no_cns_pp != 0 & p_value_no_cns_pp < 0.05,]$gene) # [1] 118

uniqueN(target_proportions[observed_prop_ap != 0,]$gene) # [1] 433
uniqueN(target_proportions[observed_prop_ap != 0 & p_value_ap < 0.05,]$gene) # [1] 129

uniqueN(target_proportions[observed_prop_qc_pp_drugbank_only != 0,]$gene) # [1] 340
uniqueN(target_proportions[observed_prop_qc_pp_drugbank_only != 0 & p_value_qc_pp_drugbank_only < 0.05,]$gene) # [1] 86

uniqueN(target_proportions[observed_prop_qc_plus_pp_drugbank_only != 0,]$gene) # [1] 256
uniqueN(target_proportions[observed_prop_qc_plus_pp_drugbank_only != 0 & p_value_qc_plus_pp_drugbank_only < 0.05,]$gene) # [1] 78

uniqueN(target_proportions[observed_prop_no_cns_pp_drugbank_only != 0,]$gene) # [1] 213
uniqueN(target_proportions[observed_prop_no_cns_pp_drugbank_only != 0 & p_value_no_cns_pp_drugbank_only < 0.05,]$gene) # [1] 44

uniqueN(target_proportions[observed_prop_ap_drugbank_only != 0,]$gene) # [1] 82
uniqueN(target_proportions[observed_prop_ap_drugbank_only != 0 & p_value_ap_drugbank_only < 0.05,]$gene) # [1] 38

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save as sheet
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
sheet4 <- copy(target_proportions)

# ==================================================================================================================================================================
#
# SHEET 5 - Overlapping targets
#
# ==================================================================================================================================================================
sheet5 <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_targets/pp_ap_overlap_significance_drugbank_only.tsv")

# ==================================================================================================================================================================
#
# SHEET 6 - Mechanism proportions
#
# ==================================================================================================================================================================

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load mechanism data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mech <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_with_actions.tsv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert mechanism data into upregular or downregulator and remove everything else
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Downregualtion
mech[action=="inhibitor" | 
		action=="antagonist" | 
		action=="blocker" | 
		action=="negative modulator" | 
		action=="inactivator" | 
		action=="suppressor" | 
		action=="weak inhibitor" | 
		action=="inhibitory allosteric modulator", action := "downregulator"]

# Upregulation
mech[action=="agonist" | 
		action=="activator" | 
		action=="inducer" | 
		action=="potentiator" | 
		action=="positive allosteric modulator" | 
		action=="positive modulator" | 
		action=="stimulator", action := "upregulator"]

# Limit to only up/downregulation
mech <- unique(mech[action == "downregulator" | action=="upregulator",])

# Add combined symbol and mechanism column
mech[,gene_mech := paste(gene,action,sep="_")]

# Set up for table
mech_proportions <- unique(mech[,.(gene_mech)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load QC Propsychotics mechs (both significant and non significant)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qc_pp_mechs_sig <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_mechanisms/propsychotic_significant_mechanisms_all.tsv")[,.(gene_mech,observed_prop_qc_pp=observed_prop,p_value_qc_pp=p_value)]
uniqueN(qc_pp_mechs_sig) # [1] 164

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load QC+ Propsychotics mechs (both significant and non significant)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qc_plus_pp_mechs_sig <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_mechanisms/propsychotic_significant_mechanisms_int.tsv")[,.(gene_mech,observed_prop_qc_plus_pp=observed_prop,p_value_qc_plus_pp=p_value)]
uniqueN(qc_plus_pp_mechs_sig) # [1] 116

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load antipsychotics mechs (both significant and non significant)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ap_mechs_sig <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/significant_mechanisms/antipsychotic_significant_mechanisms.tsv")[,.(gene_mech,observed_prop_ap=observed_prop,p_value_ap=p_value)]
uniqueN(ap_mechs_sig) # [1] 51

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mech_proportions <- merge(mech_proportions,qc_pp_mechs_sig,by="gene_mech",all.x=TRUE)
mech_proportions <- merge(mech_proportions,qc_plus_pp_mechs_sig,by="gene_mech",all.x=TRUE)
mech_proportions <- merge(mech_proportions,ap_mechs_sig,by="gene_mech",all.x=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Relabel any NA's
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sheet6 <- copy(mech_proportions)
for (i in c("observed_prop_qc_pp", "observed_prop_qc_plus_pp", "observed_prop_ap")){sheet6[is.na(get(i)), (i):=0]}
for (i in c("p_value_qc_pp", "p_value_qc_plus_pp", "p_value_ap")){sheet6[is.na(get(i)), (i):=1]}

# ==================================================================================================================================================================
#
# SHEET 7 - More readable combination of target and mechanism proportions
#
# ==================================================================================================================================================================

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simplify target data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
targets <- unique(sheet4[,.(gene,observed_prop_qc_pp_targ=observed_prop_qc_pp,observed_prop_ap_targ=observed_prop_ap,p_value_qc_pp_targ=p_value_qc_pp,p_value_ap_targ=p_value_ap)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simplify mechanism data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mechanisms <- unique(sheet6[,.(gene_mech,observed_prop_qc_pp_mech=observed_prop_qc_pp,observed_prop_ap_mech=observed_prop_ap,p_value_qc_pp_mech=p_value_qc_pp,p_value_ap_mech=p_value_ap)])
setDT(mechanisms)[, paste0("gene_mech", 1:2) := tstrsplit(gene_mech, "_")]
setnames(mechanisms,"gene_mech1", "gene")
setnames(mechanisms,"gene_mech2", "mechanism")
mechanisms[, gene_mech := NULL]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Isolate upregulation and downregulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
up <- mechanisms[mechanism == "upregulator",]
down <- mechanisms[mechanism == "downregulator",]
setnames(up,c("observed_prop_qc_pp_mech", "observed_prop_ap_mech", "p_value_qc_pp_mech", "p_value_ap_mech"), c("observed_prop_qc_pp_mech_upreg", "observed_prop_ap_mech_upreg", "p_value_qc_pp_mech_upreg", "p_value_ap_mech_upreg"))
setnames(down,c("observed_prop_qc_pp_mech", "observed_prop_ap_mech", "p_value_qc_pp_mech", "p_value_ap_mech"), c("observed_prop_qc_pp_mech_downreg", "observed_prop_ap_mech_downreg", "p_value_qc_pp_mech_downreg", "p_value_ap_mech_downreg"))
up[,mechanism := NULL]
down[,mechanism := NULL]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sheet7 <- copy(targets)
sheet7 <- merge(sheet7,up,by="gene",all.x=TRUE,all.y=TRUE)
sheet7 <- merge(sheet7,down,by="gene",all.x=TRUE,all.y=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Relabel any NA's
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in c("observed_prop_qc_pp_mech_upreg", "observed_prop_ap_mech_upreg",
			"observed_prop_qc_pp_mech_downreg", "observed_prop_ap_mech_downreg")){sheet7[is.na(get(i)), (i):=0]}

for (i in c("p_value_qc_pp_mech_upreg", "p_value_ap_mech_upreg",
			"p_value_qc_pp_mech_downreg", "p_value_ap_mech_downreg")){sheet7[is.na(get(i)), (i):=1]}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Limit to rows with at least one non-zero proportion
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sheet7 <- sheet7[observed_prop_qc_pp_mech_upreg > 0 | 
					observed_prop_ap_mech_upreg > 0 | 
					observed_prop_qc_pp_mech_downreg > 0 | 
					observed_prop_ap_mech_downreg > 0,]

setorder(sheet7,-observed_prop_qc_pp_targ)

# ==================================================================================================================================================================
#
# SHEET 8 - SCHEMA/MAGMA/PGC3 nearest gene enrichment
#
# ==================================================================================================================================================================
sheet8 <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/sz_targets/enrichment_results.tsv")

# ==================================================================================================================================================================
#
# WRITE OUT
#
# ==================================================================================================================================================================

write.xlsx(sheet1, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/scz_drug_targets_paper_tables.xlsx", sheetName="QC Propsychotic Drugs", row.names=FALSE, append=TRUE)
write.xlsx(sheet2, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/scz_drug_targets_paper_tables.xlsx", sheetName="QC+ Propsychotic Drugs", row.names=FALSE, append=TRUE)
write.xlsx(sheet3, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/scz_drug_targets_paper_tables.xlsx", sheetName="Antipsychotic Drugs", row.names=FALSE, append=TRUE)
write.xlsx(sheet4, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/scz_drug_targets_paper_tables.xlsx", sheetName="Target Proportions", row.names=FALSE, append=TRUE)
write.xlsx(sheet5, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/scz_drug_targets_paper_tables.xlsx", sheetName="Overlapping Targets", row.names=FALSE, append=TRUE)
write.xlsx(sheet6, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/scz_drug_targets_paper_tables.xlsx", sheetName="Mechanism Proportions", row.names=FALSE, append=TRUE)
write.xlsx(sheet7, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/scz_drug_targets_paper_tables.xlsx", sheetName="Targets and Mechanisms", row.names=FALSE, append=TRUE)
write.xlsx(sheet8, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/paper_tables/scz_drug_targets_paper_tables.xlsx", sheetName="SCZ Variant Enrichment", row.names=FALSE, append=TRUE)

