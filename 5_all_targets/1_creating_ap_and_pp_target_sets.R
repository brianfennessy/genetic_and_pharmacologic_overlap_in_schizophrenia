# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/targets/drugbank_seachange_targets_final.tsv")
tar[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load APs and PPs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/propsychotics_final.tsv")
ap <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/antipsychotics/ap_drugs_with_targets.tsv")

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
fwrite(pp_targs_final, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/targets/propsychotic_targets_all.tsv")
fwrite(ap_targs_final, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/targets/antipsychotic_targets_all.tsv")