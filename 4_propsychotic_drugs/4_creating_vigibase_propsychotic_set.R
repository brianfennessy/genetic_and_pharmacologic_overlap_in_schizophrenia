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
vgb <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_single_drug_and_statistical_thresholds.tsv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - RxNorm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rxn <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/rxnorm_refined.tsv")
rxn_in <- unique(rxn[tty=="IN",])
rxn_names <- unique(rxn_in[sab=="RXNORM",.(rxcui,drug_name)])
rxn_names[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load MedDRA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meddra <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/meddra_psychosis.tsv")
term_names <- unique(meddra[,.(pt_code,pt_name)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Limit VigiBase to propsychotic MedDRA terms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_pp <- vgb[pt_code %in% term_names$pt_code,]

uniqueN(vgb_pp[,.(rxcui,pt_code)]) # [1] 1025
uniqueN(vgb_pp$rxcui) # [1] 324
uniqueN(vgb_pp$pt_code) # [1] 77

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load ATC antipsychotics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ap <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/antipsychotics/ap_drugs.tsv")
ap[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exclude antipsychotics from VigiBase propsychotics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_pp <- vgb_pp[!(rxcui %in% ap$rxcui),]

uniqueN(vgb_pp[,.(rxcui,pt_code)]) # [1] 882
uniqueN(vgb_pp$rxcui) # [1] 301
uniqueN(vgb_pp$pt_code) # [1] 71

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load SIDER indications
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdr_master <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/sider_original.tsv")
sdr_ind <- sdr_master[meddra_types == "ind" & meddra %in% term_names$pt_code,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exclude SIDER indications paired with MedDRA psychotic side effects from VigiBase propsychotics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_pp_final <- unique(vgb_pp[!(rxcui %in% sdr_ind$rxcui),])

uniqueN(vgb_pp_final[,.(rxcui,pt_code)]) # [1] 732
uniqueN(vgb_pp_final$rxcui) # [1] 276
uniqueN(vgb_pp_final$pt_code) # [1] 66

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")
tar[,rxcui := as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ISOLATE - Sets of VigiBase propsychotics that have target data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vgb_pp_all <- unique(vgb_pp_final[,.(rxcui,drug_name)])
vgb_pp_w_targs <- unique(vgb_pp_final[rxcui %in% tar$rxcui,.(rxcui,drug_name)])
vgb_pp_final_w_targs <- unique(vgb_pp_final[rxcui %in% tar$rxcui,])

vgb_pp_all[drug_name=="",]
#                          rxcui drug_name
# 1:        763096|763098|763100          
# 2: 797629|797631|797633|797635          
# 3: 798262|798264|798266|798268          
# 4:               798264|798266          
# 5: 901505|901507|901509|901511   

# Note that none of the pipe-delimited rows have targets anyway so no need to stress over them, but I listed what they are below if needed to explain!
tar[rxcui=="763096" | rxcui=="763098" | rxcui=="763100" | rxcui=="797629" | rxcui=="797631" | rxcui=="797633" | rxcui=="797635" | rxcui=="798262" | rxcui=="798264" | rxcui=="798266" | rxcui=="798268" | rxcui=="798264" | rxcui=="798266" | rxcui=="901505" | rxcui=="901507" | rxcui=="901509" | rxcui=="901511",]
# Empty data.table (0 rows and 4 cols): rxcui,drug_name,gene,database

rxn_in[rxcui=="763096" & sab=="RXNORM",]
#     rxcui                                     drug_name    sab tty   code
# 1: 763096 poliovirus vaccine inactivated type 1 mahoney RXNORM  IN 763096

rxn_in[rxcui=="763098" & sab=="RXNORM",]
#     rxcui                                   drug_name    sab tty   code
# 1: 763098 poliovirus vaccine inactivated type 2 mef 1 RXNORM  IN 763098

rxn_in[rxcui=="763100" & sab=="RXNORM",]
#     rxcui                                     drug_name    sab tty   code
# 1: 763100 poliovirus vaccine inactivated type 3 saukett RXNORM  IN 763100

rxn_in[rxcui=="797629" & sab=="RXNORM",]
#     rxcui                                                                                              drug_name    sab tty   code
# 1: 797629 neisseria meningitidis serogroup a capsular polysaccharide diphtheria toxoid protein conjugate vaccine RXNORM  IN 797629

rxn_in[rxcui=="797631" & sab=="RXNORM",]
#     rxcui                                                                                              drug_name    sab tty   code
# 1: 797631 neisseria meningitidis serogroup c capsular polysaccharide diphtheria toxoid protein conjugate vaccine RXNORM  IN 797631

rxn_in[rxcui=="797633" & sab=="RXNORM",]
#     rxcui                                                                                                  drug_name    sab tty   code
# 1: 797633 neisseria meningitidis serogroup w 135 capsular polysaccharide diphtheria toxoid protein conjugate vaccine RXNORM  IN 797633

rxn_in[rxcui=="797635" & sab=="RXNORM",]
#     rxcui                                                                                              drug_name    sab tty   code
# 1: 797635 neisseria meningitidis serogroup y capsular polysaccharide diphtheria toxoid protein conjugate vaccine RXNORM  IN 797635

rxn_in[rxcui=="798262" & sab=="RXNORM",]
#     rxcui                                       drug_name    sab tty   code
# 1: 798262 l1 protein human papillomavirus type 11 vaccine RXNORM  IN 798262

rxn_in[rxcui=="798264" & sab=="RXNORM",]
#     rxcui                                       drug_name    sab tty   code
# 1: 798264 l1 protein human papillomavirus type 16 vaccine RXNORM  IN 798264

rxn_in[rxcui=="798266" & sab=="RXNORM",]
#     rxcui                                       drug_name    sab tty   code
# 1: 798266 l1 protein human papillomavirus type 18 vaccine RXNORM  IN 798266

rxn_in[rxcui=="798268" & sab=="RXNORM",]
#     rxcui                                      drug_name    sab tty   code
# 1: 798268 l1 protein human papillomavirus type 6 vaccine RXNORM  IN 798268

rxn_in[rxcui=="798264" & sab=="RXNORM",]
#     rxcui                                       drug_name    sab tty   code
# 1: 798264 l1 protein human papillomavirus type 16 vaccine RXNORM  IN 798264

rxn_in[rxcui=="798266" & sab=="RXNORM",]
#     rxcui                                       drug_name    sab tty   code
# 1: 798266 l1 protein human papillomavirus type 18 vaccine RXNORM  IN 798266

rxn_in[rxcui=="901505" & sab=="RXNORM",]
#     rxcui                                                                                      drug_name    sab tty   code
# 1: 901505 neisseria meningitidis serogroup a oligosaccharide diphtheria crm197 protein conjugate vaccine RXNORM  IN 901505

rxn_in[rxcui=="901507" & sab=="RXNORM",]
#     rxcui                                                                                      drug_name    sab tty   code
# 1: 901507 neisseria meningitidis serogroup c oligosaccharide diphtheria crm197 protein conjugate vaccine RXNORM  IN 901507

rxn_in[rxcui=="901509" & sab=="RXNORM",]
#     rxcui                                                                                          drug_name    sab tty   code
# 1: 901509 neisseria meningitidis serogroup w 135 oligosaccharide diphtheria crm197 protein conjugate vaccine RXNORM  IN 901509

rxn_in[rxcui=="901511" & sab=="RXNORM",]
#     rxcui                                                                                      drug_name    sab tty   code
# 1: 901511 neisseria meningitidis serogroup y oligosaccharide diphtheria crm197 protein conjugate vaccine RXNORM  IN 901511

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(vgb_pp_final, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_propsychotics_all_information.tsv")
fwrite(vgb_pp_final_w_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_propsychotics_all_information_with_targets.tsv")
fwrite(vgb_pp_all, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_propsychotics_drugs_only.tsv")
fwrite(vgb_pp_w_targs, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/vigibase_propsychotics_drugs_with_targets_only.tsv")


