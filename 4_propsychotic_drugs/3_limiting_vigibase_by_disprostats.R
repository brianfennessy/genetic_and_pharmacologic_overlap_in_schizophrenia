# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load VigiBase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_filter.tsv")

uniqueN(vgb[,.(rxcui,pt_code)]) # [1] 1024552
uniqueN(vgb$rxcui) # [1] 3721
uniqueN(vgb$pt_code) # [1] 16232

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Limit VigiBase to prr >= 3, nreport.with.drug.with.ae >= 3 and chisq.yates.stat >= 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_dispro <- vgb[prr >= 3 & chisq.yates.stat >= 4 & nreport.with.drug.with.ae >= 3,]

uniqueN(vgb_dispro[,.(rxcui,pt_code)]) # [1] 81675
uniqueN(vgb_dispro$rxcui) # [1] 2667
uniqueN(vgb_dispro$pt_code) # [1] 9035

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(vgb_dispro, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_and_statistical_thresholds.tsv")
