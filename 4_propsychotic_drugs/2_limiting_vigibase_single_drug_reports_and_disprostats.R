# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Store user home directory and navigate to it
user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
setwd(user_home_directory)

# Get the official home directory of the user
home_dir <- getwd()

# Navigate to the raw data zip file and unzip
setwd(home_dir)

source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load VigiBase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mydatout <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_original.tsv"
mydatstat <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_original_report_stats.tsv"
mydat <- fread(mydatout)
mystats <- fread(mydatstat)
uniqueN(mydat) # [1] 148056551

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define strict set, get stats for number of drugs per report, number of drug-ae pairs per report
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strict.reports <- mystats[n.ae==1 & n.ing==1 & n.ing.mapped==1 & n.ind.mapped==1]$UMCReportId
mydat.strict.drugae <- unique(mydat[UMCReportId %in% strict.reports,.(drug=rxcui_output, ae=adr.meddra, strict=1)])
report.counts <- mystats[,.(report=UMCReportId,report_drug_count=n.ing,report_drug_ae_pairs_count=n.ing.ae.pairs)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# DISPROSTATS - FILTERED VIGIBASE - single drug reports only
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up input for disprostats and run disprostats
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
badreports <- c()
badreports <- unique(c( badreports, mystats[n.ing!=1]$UMCReportId))
dp.filt.inp <- unique(mydat[!UMCReportId %in% badreports,.(report=UMCReportId, drug=rxcui_output, ae=adr.meddra)][!is.na(drug) & !is.na(ae) & !is.na(report)])
uniqueN(dp.filt.inp) # [1] 13779455
dp.filt.out <- disprostats(dp.filt.inp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save a file in case something gets messed up down the line
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dp.filt.out.save <- copy(dp.filt.out)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Resume merging
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uniqueN(dp.filt.out) # [1] 60399272
dp.filt.out <- merge( unique(dp.filt.inp[,.(drug,ae,report)]), dp.filt.out, by=c("drug","ae") )
uniqueN(dp.filt.out) # [1] 13779455

dp.filt.out <- merge(dp.filt.out, mydat.strict.drugae, all.x=TRUE)
uniqueN(dp.filt.out) # [1] 13779455

dp.filt.out[is.na(strict),strict:=0]
dp.filt.out[, propreport.with.drug.with.ae := nreport.with.drug.with.ae/(nreport.with.drug.with.ae+nreport.with.drug.without.ae) ]
dp.filt.out[, chisq.yates.p.adj.fdr := p.adjust(chisq.yates.p, method="fdr") ]
dp.filt.out[, chisq.yates.p.adj.bon := p.adjust(chisq.yates.p, method="bonferroni") ]

uniqueN(dp.filt.out) # [1] 13779455

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Backup for merging drug and side effect information
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_sd_with_reports <- unique(dp.filt.out[,.(report,rxcui=drug,pt_code=ae,nreport.with.drug.with.ae,nreport.with.drug.without.ae,
							nreport.with.ae.without.drug,nreport.without.drug.without.ae,ror,ror_se,ror_95ci_ul,
							ror_95ci_ll,prr,prr_se,prr_95ci_ul,prr_95ci_ll,IC,IC_var,IC_sd,chisq.yates.stat,chisq.yates.p,strict,
							propreport.with.drug.with.ae,chisq.yates.p.adj.fdr,chisq.yates.p.adj.bon)])

vgb_sd_no_reports <- unique(dp.filt.out[,.(rxcui=drug,pt_code=ae,nreport.with.drug.with.ae,nreport.with.drug.without.ae,
							nreport.with.ae.without.drug,nreport.without.drug.without.ae,ror,ror_se,ror_95ci_ul,
							ror_95ci_ll,prr,prr_se,prr_95ci_ul,prr_95ci_ll,IC,IC_var,IC_sd,chisq.yates.stat,chisq.yates.p,strict,
							propreport.with.drug.with.ae,chisq.yates.p.adj.fdr,chisq.yates.p.adj.bon)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load RxNorm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rxn <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/rxnorm_refined.tsv")
rxn_in <- unique(rxn[tty=="IN",])
rxn_names <- unique(rxn_in[sab=="RXNORM",.(rxcui,drug_name)])
rxn_names[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load MedDRA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meddra <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/meddra_original.tsv")
term_names <- unique(meddra[,.(pt_code,pt_name)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ADD - MedDRA and RxNorm information to VigiBase propsychotics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_sd_no_reports <- merge(vgb_sd_no_reports,term_names,by="pt_code",all.x=TRUE)
vgb_sd_no_reports <- merge(vgb_sd_no_reports,rxn_names,by="rxcui",all.x=TRUE)

vgb_sd_with_reports <- merge(vgb_sd_with_reports,term_names,by="pt_code",all.x=TRUE)
vgb_sd_with_reports <- merge(vgb_sd_with_reports,rxn_names,by="rxcui",all.x=TRUE)

vgb_sd_no_reports <- unique(vgb_sd_no_reports[,.(rxcui,drug_name,pt_code,pt_name,nreport.with.drug.with.ae,nreport.with.drug.without.ae,
							nreport.with.ae.without.drug,nreport.without.drug.without.ae,ror,ror_se,ror_95ci_ul,
							ror_95ci_ll,prr,prr_se,prr_95ci_ul,prr_95ci_ll,IC,IC_var,IC_sd,chisq.yates.stat,chisq.yates.p,strict,
							propreport.with.drug.with.ae,chisq.yates.p.adj.fdr,chisq.yates.p.adj.bon)])

vgb_sd_with_reports <- unique(vgb_sd_with_reports[,.(report,rxcui,drug_name,pt_code,pt_name,nreport.with.drug.with.ae,nreport.with.drug.without.ae,
							nreport.with.ae.without.drug,nreport.without.drug.without.ae,ror,ror_se,ror_95ci_ul,
							ror_95ci_ll,prr,prr_se,prr_95ci_ul,prr_95ci_ll,IC,IC_var,IC_sd,chisq.yates.stat,chisq.yates.p,strict,
							propreport.with.drug.with.ae,chisq.yates.p.adj.fdr,chisq.yates.p.adj.bon)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")
tar[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE - Separate lists of VigiBase drugs with and without targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_sd_report_drugs_all <- unique(vgb_sd_with_reports[,.(rxcui,drug_name)])
vgb_sd_report_drugs_w_targs <- unique(vgb_sd_with_reports[rxcui %in% tar$rxcui,.(rxcui,drug_name)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE - out VigiBase limited to single drug reports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(vgb_sd_with_reports, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_filter_with_report_number.tsv")
fwrite(vgb_sd_no_reports, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_filter.tsv")
fwrite(vgb_sd_report_drugs_all, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_filter_drugs_only.tsv")
fwrite(vgb_sd_report_drugs_w_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_filter_with_targets_drugs_only.tsv")



