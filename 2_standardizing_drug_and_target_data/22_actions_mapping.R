# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
library(tidyverse)
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load DB gene ID to gene name
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene_map <- fread("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/dbgeneid2uni.tsv", sep="\t")
gene_map[,gene_other:=NULL]
gene_map <- gene_map[!is.na(gene),]
gene_map <- unique(gene_map)

uniqueN(gene_map$uniprot) # [1] 4378
uniqueN(gene_map$db_gene) # [1] 4373
uniqueN(gene_map$gene) # [1] 3981
uniqueN(gene_map[,.(uniprot,db_gene)]) # [1] 4597
uniqueN(gene_map[,.(db_gene,gene)]) # [1] 4595
uniqueN(gene_map[,.(uniprot,gene)]) # [1] 4378
uniqueN(gene_map[,.(uniprot,db_gene,gene)]) # [1] 4597

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load DB drug ID to rxn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
drug_map <- fread("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/dbdrugid2rxn.tsv", sep="\t")

uniqueN(drug_map$dbid_drug) # [1] 2150
uniqueN(drug_map$rxcui) # [1] 2141
uniqueN(drug_map[,.(dbid_drug,rxcui)]) # [1] 2196

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load rxnorm ingredients
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ing_map <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_study_2020/rxnorm_files/rxn_ingredient_mapping.tsv", sep="\t")
ing_map <- unique(ing_map[,.(rxcui,drug_name=drug)])
ing_map[,rxcui:=as.character(rxcui)]

uniqueN(ing_map$rxcui) # [1] 11967
uniqueN(ing_map$drug_name) # [1] 11961
uniqueN(ing_map[,.(rxcui,drug_name)]) # [1] 11967

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add rxnorm ingredient names to drug mapping
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
drug_map <- merge(drug_map,ing_map,by="rxcui",all.x=TRUE)

# Ask Alex and Doug what to do with these, for now I'll keep them in
drug_map[is.na(drug_name),]
#          rxcui dbid_drug drug_name
# 1:  10995|4582   DB09327      <NA> # Uracil | Tegafur
# 2: 227518|6383   DB00032      <NA> # Luteinizing Hormone | Follicle Stimulating Hormone

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load actions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
actions_path <- "/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_ACTIONS.rrf"
actions <- fread(actions_path)[,.(dbid_drug=V1,db_gene=V3,action=V4)]

uniqueN(actions$dbid_drug) # [1] 7399
uniqueN(actions$db_gene) # [1] 4865
uniqueN(actions[,.(dbid_drug, db_gene)]) # [1] 18495
uniqueN(actions) # [1] 18495

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge like it's your middle name
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add genetic information
final <- merge(actions,gene_map,by="db_gene",all.x=TRUE)

# Add drug information
final <- merge(final,drug_map,by="dbid_drug",all.x=TRUE)

# Get rid of NA's
final <- final %>% drop_na()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create final mapping file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final <- unique(final[,.(rxcui,drug_name,gene,action)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How'd we do? --> Keep in mind that this isn't for our targets, this is just for the actions table in general!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

uniqueN(actions[,.(dbid_drug,db_gene,action)]) # [1] 18495
uniqueN(final[,.(rxcui,gene,action)]) # [1] 10563

uniqueN(final$rxcui) # [1] 2122
uniqueN(final$drug_name) # [1] 2122
uniqueN(final$gene) # [1] 2760
uniqueN(final[,.(rxcui,drug_name,gene)]) # [1] 10359

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(final, sep="\t", "/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/mapped_target_actions.tsv")