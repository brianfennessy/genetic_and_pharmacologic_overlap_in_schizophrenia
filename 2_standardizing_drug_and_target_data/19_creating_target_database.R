# ============================================================================================================================================================================
# ============================================================================================================================================================================
#
# DRUGBANK
#
# ============================================================================================================================================================================
# ============================================================================================================================================================================

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Drug id maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_rxn2rxn.Rd")
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_db2rxn_2020.Rd")
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_sc2rxn.Rd")
db2rxn$db.mapped[,rxcui:=as.character(rxcui)]
setnames(db2rxn$db, "DrugBank ID", "DBID")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Drugbank target data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
targets_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/target-all-uniprot-links.csv"
dbk <- fread(targets_path)[,.(dbid=`DrugBank ID`,uniprot=`UniProt ID`, type="target")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - RxNorm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/rxnorm_original.Rdata")
rxnLatest.conso[, STR.mod := makemod(STR) ]
rxn_refined <- unique(rxnLatest.conso[,.(rxcui=RXCUI,drug_name=STR.mod,sab=SAB,tty=TTY,code=CODE)])
fwrite(rxn_refined, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/rxnorm_refined.tsv")

rxn <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/propsychotics/rxnorm_refined.tsv", sep="\t")
rxn[,rxcui:=as.character(rxcui)]
rxn_in <- rxn[tty=="IN" & sab=="RXNORM",.(rxcui,drug_name)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here, we isolate the major drugbank file to the drugbank ID, whether or not that ID is mapped to rxn, and whether or not that ID is approved
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp0 <- db2rxn$db[,.(dbid=DBID, rxnmapped.rxn, status.approved)] # all drugbank drugs and their status in drugbank and if we mapped them to rxn

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here, we get...
# A unque list of drugbank-rxnorm mapped id codes merged with an rxnorm ingredient mapping file (rxn2rxn)
# And then isolates the drug bank id's that are not NAs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dbid.with.final.rxcui <- unique(merge(db2rxn$db.mapped[,.(dbid, rxcui_input=rxcui)], rxn2rxn[[1]])[!is.na(rxcui_output)]$dbid)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here, we change columns in tmp0 based on whether or not the drug bank ID is mapped to rxnorm or not
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp0[, rxnmapped.rxn:="not.mapped.to.rxn.sab"]
tmp0[dbid %in% dbid.with.final.rxcui, rxnmapped.rxn:="mapped.to.rxn.sab"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here, we make a new data table of all drug bank ID's and start with them all "having targets"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp1 <- unique(dbk[,.(dbid, has_targets=1)]) #drugbank drugs with targets

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Then we take tmp 1 and ID all the drug bank ID's that do not have targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp2 <- merge(tmp0, tmp1, all=T) 
tmp2[is.na(has_targets), has_targets:=0]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here, we...
# 1) Read in the mapping file for uniprot to gene symbol
# 2) Create a drugbank table with drug bank ID and uniprot pairs mapped to rxcui's
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dbkmap <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drugbank_uniprot_id_map.tsv", na="", col.names = c("uniprot","gene","gene_other")) 
dbk_mapping <- merge(dbk, merge(db2rxn$db.mapped[,.(dbid, rxcui_input=rxcui)], rxn2rxn[[1]]), by="dbid")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here, we complete the mapping process
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dbk_final <- unique(merge(dbk_mapping, dbkmap, by="uniprot")[,.(rxcui=rxcui_output, gene, gene_other, database="drugbank")])
# dbk_final <- unique(merge(dbk_mapping, dbkmap, by="uniprot")
	
# ============================================================================================================================================================================
# ============================================================================================================================================================================
#
# SEACHANGE
#
# ============================================================================================================================================================================
# ============================================================================================================================================================================

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD - SeaChange files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/seachange")
scd0 <- fread("SeaChange_SC_binding_cfp_1-13-14.csv", header=T,sep=",", na.strings="NULL")
scd1 <- fread("SeaChange_SC_binding_path_1-13-14.csv", header=T,sep=",", na.strings="NULL")
scd2 <- fread("SeaChange_SC_functional_cfp_1-13-14.csv", header=T,sep=",", na.strings="NULL")
scd3 <- fread("SeaChange_SC_functional_path_1-13-14.csv",header=T,sep=",", na.strings="NULL")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY - the SeaChange method
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scd0 <- scd0[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "binding_cfp"] 
scd1 <- scd1[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "binding_path"]
scd2 <- scd2[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "functional_cfp"] 
scd3 <- scd3[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "functional_path"] 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINE - all SeaChange files and reformat them
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sea <- rbind( scd0,scd1,scd2,scd3,use.names=T)
sea <- unique(sea[,.(chembl_id=ChEMBL, initial_source="seachange", uniprot_title=Mol, ncbi_gene_id=GeneID,
                      seachange_pval=Pvalue, seachange_maxtc=Max_Tc, seachange_method)][!is.na(ncbi_gene_id),])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ORGANIZE - columns
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fwrite(unique(sea[,.(uniprot_title)]), row=F, quo=F, col=F, file="~/www/files/sc_up_ids") 
# ##
# ##map at https://www.uniprot.org/uploadlists/ then read in results
# ##
cnames <- c("uniprot","uniprot_title","uniprot_species","gene","gene_other")
seamap <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/seachange_uniprot_id_map.tsv", na="", col.names = cnames)
sea <- merge(sea, seamap, by="uniprot_title")
sea <- unique(sea[,.(chembl_id, gene, gene_other)])
tmp <- merge( sc2rxn[,.(chembl_id,rxcui_input=as.character(rxcui))], rxn2rxn[[1]], by="rxcui_input")[,.(chembl_id,rxcui=rxcui_output)]
sea <- unique(merge(sea, tmp)[,.(rxcui, gene, gene_other, database="seachange")])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINE - drugbank and seachange
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar <- rbind(dbk_final, sea)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CURATE - final data set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tar[,gene:=toupper(gene)]
tar[,gene_other:=toupper(gene_other)]
tar[,rxcui:=as.integer(rxcui)] #removes multiingred drugs
tar <- tar[!is.na(gene),] # There are some NA values for genes
tar <- tar[!is.na(rxcui),] # There are some NA values for drugs
tar <- unique(tar)
tar[,rxcui:=as.character(rxcui)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ADD - RxNorm Ingredient Names
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final <- merge(tar,rxn_in,by="rxcui",all.x=TRUE)[,.(rxcui,drug_name,gene,gene_other,database)]
final <- final[gene != "NA",]

final_w_gene_other <- unique(final)
final_simplified <- unique(final[,.(rxcui,drug_name,gene,database)])

# ============================================================================================================================================================================
# ============================================================================================================================================================================
#
# WRITE OUT
#
# ============================================================================================================================================================================
# ============================================================================================================================================================================

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE - final data set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(final_w_gene_other, row=F, quo=F, sep="\t", file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_with_gene_other.tsv")
fwrite(final_simplified, row=F, quo=F, sep="\t", file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/mapped_files/targets/drugbank_seachange_targets_final.tsv")


