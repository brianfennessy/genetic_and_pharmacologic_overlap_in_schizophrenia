# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
options(width=160)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/target-all-uniprot-links.csv"
dbk <- fread(target_path)[,.(dbid=`DrugBank ID`,uniprot=`UniProt ID`, type="target")]

uniqueN(dbk$dbid) # [1] 7269
uniqueN(dbk$uniprot) # [1] 4766
uniqueN(dbk[,.(dbid, uniprot)]) # [1] 19317

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load actions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
actions_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_ACTIONS.rrf"
actions <- fread(actions_path)[,.(dbid=V1,dbid_gene=V3,action=V4)]

uniqueN(actions$dbid) # [1] 7399
uniqueN(actions$dbid_gene) # [1] 4865
uniqueN(actions[,.(dbid, dbid_gene)]) # [1] 18495

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Alex's mapping code from drug2target
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Drug id maps
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_rxn2rxn.Rd")
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_db2rxn_2020.Rd")
load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_sc2rxn.Rd")
db2rxn$db.mapped[,rxcui:=as.character(rxcui)]

# Here, we isolate the major drugbank file to the drugbank ID, whether or not that ID is mapped to rxn, and whether or not that ID is approved
tmp0 = db2rxn$db[,.(dbid=`DrugBank ID`, rxnmapped.rxn, status.approved)] # all drugbank drugs and their status in drugbank and if we mapped them to rxn

# Here, we get...
# A unque list of drugbank-rxnorm mapped id codes merged with an rxnorm ingredient mapping file (rxn2rxn)
# And then isolates the drug bank id's that are not NAs
dbid.with.final.rxcui = unique(merge(db2rxn$db.mapped[,.(dbid, rxcui_input=rxcui)], rxn2rxn[[1]])[!is.na(rxcui_output)]$dbid)

# Here, we change columns in tmp0 based on whether or not the drug bank ID is mapped to rxnorm or not
tmp0[, rxnmapped.rxn:="not.mapped.to.rxn.sab"]
tmp0[dbid %in% dbid.with.final.rxcui, rxnmapped.rxn:="mapped.to.rxn.sab"]

# Here, we make a new data table of all drug bank ID's and start with them all "having targets"
tmp1 = unique(dbk[,.(dbid, has_targets=1)]) # drugbank drugs with targets

# Then we take tmp 1 and ID all the drug bank ID's that do not have targets
tmp2 = merge(tmp0, tmp1, all=T) 
tmp2[is.na(has_targets), has_targets:=0]

# Here, we...
# 1) Read in the mapping file for uniprot to gene symbol
# 2) Create a drugbank table with drug bank ID and uniprot pairs mapped to rxcui's
dbkmap=fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drugbank_uniprot_id_map.tsv", na="", col.names = c("uniprot","gene","gene_other")) 
dbk_mapping = merge(dbk, merge(db2rxn$db.mapped[,.(dbid, rxcui_input=rxcui)], rxn2rxn[[1]]), by="dbid")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create final mapping file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final <- unique(dbk_mapping[,.(dbid_drug=dbid,rxcui=rxcui_output)])

uniqueN(final$dbid_drug) # [1] 2150
uniqueN(final$rxcui) # [1] 2141
uniqueN(final[,.(dbid_drug,rxcui)]) # [1] 2196

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(final, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/dbdrugid2rxn.tsv")
