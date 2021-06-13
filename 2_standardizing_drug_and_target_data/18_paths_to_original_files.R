# Load libraries
source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")

# VigiBase - Original
cp /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv /sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_original.tsv
cp /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat_report_stats.tsv /sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_original_report_stats.tsv

mydatout <- "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_original.tsv" # Originally called "mydat"
mydatstat <- "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/vigibase_original_report_stats.tsv" # Originally called "mystats"

mydat <- fread(mydatout)
mystats <- fread(mydatstat)

# SIDER - Original
cp /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/mysdr_final_map.tsv /sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/sider_original.tsv
sdrf <- "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/sider_original.tsv"
mysdr <- fread(sdrf)

# ATC - Original
atc_path <- "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/antipsychotics/atc_original.tsv"
atc <- fread(atc_path)

# MedDRA - Original
meddra_path <- "/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/meddra_original.tsv"
meddra <- fread("meddra_path")

# RxNorm - First two lines are the same files as the bottom two, just copied
# load("/sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata")
# rxnLatest.conso[, STR.mod := makemod(STR) ]
rxn_original <- load("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/rxnorm_original.Rdata")
rxn_refined <- fread("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/propsychotics/rxnorm_refined.tsv", sep="\t")

# DrugBank - Original
targets_path = "/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/target-all-uniprot-links.csv"
dbk = fread(targets_path)[,.(dbid=`DrugBank ID`,uniprot=`UniProt ID`, type="target")]

# SeaChange - Original
setwd("/sc/arion/projects/psychgen/methods/rx/data/drug_target/seachange")
scd0 = fread("SeaChange_SC_binding_cfp_1-13-14.csv", header=T,sep=",", na.strings="NULL")
scd1 = fread("SeaChange_SC_binding_path_1-13-14.csv", header=T,sep=",", na.strings="NULL")
scd2 = fread("SeaChange_SC_functional_cfp_1-13-14.csv", header=T,sep=",", na.strings="NULL")
scd3 = fread("SeaChange_SC_functional_path_1-13-14.csv",header=T,sep=",", na.strings="NULL")
scd0 = scd0[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "binding_cfp"] 
scd1 = scd1[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "binding_path"]
scd2 = scd2[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "functional_cfp"] 
scd3 = scd3[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "functional_path"] 
sea = rbind( scd0,scd1,scd2,scd3,use.names=T)
sea = unique(sea[,.(chembl_id=ChEMBL, initial_source="seachange", uniprot_title=Mol, ncbi_gene_id=GeneID,
                      seachange_pval=Pvalue, seachange_maxtc=Max_Tc, seachange_method)][!is.na(ncbi_gene_id),])