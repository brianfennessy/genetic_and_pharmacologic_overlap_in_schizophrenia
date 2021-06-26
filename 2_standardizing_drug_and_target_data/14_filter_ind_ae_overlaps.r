# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/mysdr_final_map.tsv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/mysdr_final_map.tsv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load sider
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdrf="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/mysdr_final_map.tsv"
  mysdr = fread(sdrf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# remove instances where traits is both an AE and IND for same drug
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mysdr.collapsed = mysdr[, list(meddra_types = paste(unique(meddra_type), collapse="|")), by=list(rxcui,meddra,meddra_level)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  write.table(mysdr.collapsed, row=F, quo=F, sep="\t", file=sdrf)
