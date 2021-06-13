# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
library(tidyverse)
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load targets table
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
targets_master <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/targets/drugbank_seachange_targets_final.tsv", sep="\t")
targets_master[,rxcui:=as.character(rxcui)]
tar <- targets_master[database=="drugbank"]

uniqueN(tar$rxcui) # [1] 2122
uniqueN(tar$drug_name) # [1] 2122
uniqueN(tar$gene) # [1] 2760
uniqueN(tar[,.(rxcui,drug_name,gene)]) # [1] 10359

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load mapped actions table
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
actions <- fread("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/mapped_target_actions.tsv", sep="\t")

uniqueN(actions$rxcui) # [1] 2122
uniqueN(actions$drug_name) # [1] 2122
uniqueN(actions$gene) # [1] 2760
uniqueN(actions[,.(rxcui,drug_name,gene)]) # [1] 10359

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add merging columns
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar[,to_merge := paste(rxcui,gene,sep=",")]
actions[,to_merge := paste(rxcui,gene,sep=",")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge like it's your middle name
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final <- merge(tar,actions,by="to_merge",all.x=TRUE)
final <- unique(final[,.(rxcui=rxcui.x,drug_name=drug_name.x,gene=gene.x,action,database)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How'd we do?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uniqueN(tar[,.(rxcui,gene)]) # [1] 10359
uniqueN(final[,.(rxcui,gene)]) # [1] 10359

uniqueN(final$rxcui) # [1] 2122
uniqueN(final$drug_name) # [1] 2122
uniqueN(final$gene) # [1] 2760
uniqueN(final[,.(rxcui,drug_name,gene)]) # [1] 10359

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(final, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/targets/drugbank_seachange_targets_with_actions.tsv")

