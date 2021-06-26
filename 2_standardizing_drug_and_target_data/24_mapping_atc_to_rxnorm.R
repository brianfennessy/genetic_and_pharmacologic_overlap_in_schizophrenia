# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load ATC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atc <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/atc_original.tsv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load RxNorm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rxn_master <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/rxnorm_refined.tsv")
rxn <- rxn_master[tty=="IN" & sab=="ATC",]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atc_mapped <- merge(atc,rxn,by.x="level5",by.y="code",all.x=TRUE)
unique(atc_mapped[is.na(rxcui),]$who_name)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(atc_mapped, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/atc_rxnorm_mapped.tsv")
