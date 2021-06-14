# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load VigiBase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_single_drug_filter.tsv")

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
fwrite(vgb_dispro, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_single_drug_and_statistical_thresholds.tsv")