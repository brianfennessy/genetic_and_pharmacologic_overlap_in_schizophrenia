# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load VigiBase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_propsychotics_drugs_only.tsv")
vgb_drugs <- unique(vgb$rxcui)

vgb_w_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_propsychotics_drugs_with_targets_only.tsv")
vgb_drugs_w_targs <- unique(vgb_w_targs$rxcui)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load SIDER
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdr <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/sider_propsychotics.tsv")
sdr_drugs <- unique(sdr$rxcui)

sdr_w_targs <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/sider_propsychotics_with_targets.tsv")
sdr_drugs_w_targs <- unique(sdr_w_targs$rxcui)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Isolate intersect between two sets (with targets)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
propsychotics_w_targs <- intersect(vgb_drugs_w_targs, sdr_drugs_w_targs)
pp_w_targs <- as.data.table(propsychotics_w_targs)
pp_w_targs <- pp_w_targs[,.(rxcui=propsychotics_w_targs)]
pp_w_targs[,rxcui := as.character(rxcui)]

# Check intersect in general
length(intersect(vgb_drugs,sdr_drugs)) # [1] 148
length(intersect(vgb_drugs_w_targs,sdr_drugs_w_targs)) # [1] 147

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load RxNorm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rxn <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/rxnorm_refined.tsv")
rxn_in <- unique(rxn[tty=="IN",])
rxn_names <- unique(rxn_in[sab=="RXNORM",.(rxcui,drug_name)])
rxn_names[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add drug names
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pp <- merge(pp_w_targs,rxn_names,by="rxcui",all.x=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(pp, sep="\t", "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/propsychotics_final.tsv")
