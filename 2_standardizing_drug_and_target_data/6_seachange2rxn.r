# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Seachange was delivered to us with drugs identified as chembl IDs.
# Mapping chembl to RxNorm requires intermediate IDs (i.e., chembl is not an RxNorm vocabulary). 
# Note: CHEBI IDs in seachange file to not appear correct - DONT USE THEM! 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_db2rxn.Rd
# 2. /sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata
# 3. /sc/arion/projects/psychgen/methods/rx/data/drug_target/seachange/SeaChange_SC_binding_cfp_1-13-14.csv
# 4. /sc/arion/projects/psychgen/methods/rx/data/drug_target/seachange/SeaChange_SC_binding_path_1-13-14.csv
# 5. /sc/arion/projects/psychgen/methods/rx/data/drug_target/seachange/SeaChange_SC_functional_cfp_1-13-14.csv
# 6. /sc/arion/projects/psychgen/methods/rx/data/drug_target/seachange/SeaChange_SC_functional_path_1-13-14.csv
# 7. /sc/arion/projects/psychgen/methods/rx/data/chembl/src1src22.txt
# 8. /sc/arion/projects/psychgen/methods/rx/data/chembl/src1src2.txt
# 9. /sc/arion/projects/psychgen/methods/rx/data/chembl/src1src14.txt
# 10. /sc/arion/projects/psychgen/methods/rx/data/uniII/UNIIs_10Nov2016_Records.txt
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_sc2rxn.Rd
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
  library(readr) 
  library(stringr)
  setwd("/sc/arion/projects/psychgen/methods/rx/data/drug_target/seachange")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in drugbank (its helpful here)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("/sc/arion/projects/psychgen/methods/rx/files/drug_id_map_db2rxn.Rd")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in latest rxnorm (9/4/18 release)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("/sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata")
  rxnLatest.conso[, STR.mod := makemod(STR) ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in seachange and format
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scd0 = fread("SeaChange_SC_binding_cfp_1-13-14.csv",    header=T,sep=",", na.strings="NULL")[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "binding_cfp"] 
  scd1 = fread("SeaChange_SC_binding_path_1-13-14.csv",   header=T,sep=",", na.strings="NULL")[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "binding_path"]
  scd2 = fread("SeaChange_SC_functional_cfp_1-13-14.csv", header=T,sep=",", na.strings="NULL")[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "functional_cfp"] 
  scd3 = fread("SeaChange_SC_functional_path_1-13-14.csv",header=T,sep=",", na.strings="NULL")[ , Pvalue := 10^log10_Pvalue ][, seachange_method := "functional_path"] 
  sea = rbind( scd0,scd1,scd2,scd3,use.names=T)
  sea = unique(sea[,.(chembl_id=ChEMBL, initial_source="seachange", uniprot_title=Mol, ncbi_gene_id=GeneID,
                              seachange_pval=Pvalue, seachange_maxtc=Max_Tc, seachange_method)][!is.na(ncbi_gene_id),])
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# map chembl to other ids via unichem (unii, pubchem, etc; for sources of these mapping files,see /sc/arion/projects/psychgen/methods/rx/data/chembl/00README)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  chembl2pc = fread('/sc/arion/projects/psychgen/methods/rx/data/chembl/src1src22.txt', header=T, col.names=c("chembl_id","intermediate_id"))
  chembl2db = fread('/sc/arion/projects/psychgen/methods/rx/data/chembl/src1src2.txt' , header=T, col.names=c("chembl_id","intermediate_id"))
  chembl2unii = fread('/sc/arion/projects/psychgen/methods/rx/data/chembl/src1src14.txt', header=T, col.names=c("chembl_id","intermediate_id"))
  chembl2name = fread("chembl_genericName_dump.txt", header=T)[,.(intermediate_id=makemod(pref_name),chembl_id)] #seachange file linking names to chembl
  chembl2pc = chembl2pc[chembl_id %in% unique(sea$chembl_id),.(chembl_id,intermediate_id,intermediate_source="pubchem")]
  chembl2db = chembl2db[chembl_id %in% unique(sea$chembl_id),.(chembl_id,intermediate_id,intermediate_source="drugbank")]
  chembl2unii = chembl2unii[chembl_id %in% unique(sea$chembl_id),.(chembl_id,intermediate_id,intermediate_source="unii")]
  chembl2name = chembl2name[chembl_id %in% unique(sea$chembl_id),.(chembl_id,intermediate_id,intermediate_source="seachange_pref_name")]
  mysea = rbind( chembl2pc, chembl2db, chembl2unii, chembl2name )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of drugs in seachange
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uniqueN(sea$chembl_id) #[1] 1058 ... number of drugs in seachange ... confirmed 13apr2021
  uniqueN(mysea$chembl_id) #[1] 1039 ... number mapped to an intermediate id (pubchem/drugbank/unii/chembl_name) ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# track mappings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mysea.mapped = c()
  mysea.unmapped = c()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# exact string matching of seachange names to rxnorm names (with sab "rxnorm")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scnames = makemod(unique(mysea[intermediate_source=="seachange_pref_name"]$intermediate_id))
  x1 = unique(rxnLatest.conso[STR.mod %in% scnames]$RXCUI)
  x2 = unique(rxnLatest.conso[ RXCUI %in% x1 & SAB=="RXNORM"]$RXCUI)
  x3 = unique(rxnLatest.conso[RXCUI %in% x2 & STR.mod %in% scnames,.(RXCUI,intermediate_id=STR.mod)])
  x4 = merge(x3, mysea)
  newly.mapped = x4[,.(chembl_id, rxcui=RXCUI)]
  mysea.mapped = rbind(mysea.mapped, newly.mapped)[,.(chembl_id,rxcui)]
  mysea.unmapped = rbind(mysea.unmapped, mysea[! chembl_id %in% unique(mysea.mapped$chembl_id)])
  uniqueN(mysea.mapped$chembl_id) #[1] 776 ... confirmed 13apr2021
  uniqueN(mysea.unmapped$chembl_id) #[1] 263 ... sanity check 776+263=1039... ok, good  ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sc to rxn via unii - file downloaded from unii
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  unii = fread('/sc/arion/projects/psychgen/methods/rx/data/uniII/UNIIs_10Nov2016_Records.txt', sep='\t', na.strings=c("NA","") )
  unii$PT = tolower(unii$PT)
  scmap0 = merge(mysea.unmapped,unii, by.x="intermediate_id", by.y="PT")[,.(chembl_id,intermediate_id,intermediate_source,rxcui=RXCUI)] #map by name
  scmap1 = merge(mysea.unmapped,unii, by.x="intermediate_id", by.y="UNII")[,.(chembl_id,intermediate_id,intermediate_source,rxcui=RXCUI)] #map by unii
  newly.mapped = unique(rbind(scmap0,scmap1,use.names=T))[!is.na(rxcui),.(chembl_id,rxcui)]
  mysea.mapped = rbind(mysea.mapped, newly.mapped )
  mysea.unmapped = mysea[! chembl_id %in% unique(mysea.mapped$chembl_id)]
  uniqueN(mysea.mapped$chembl_id) #[1] 794 ... confirmed 13apr2021
  uniqueN(mysea.unmapped$chembl_id) #[1] 245 ... confirmed 13apr2021
  sc2rxn = unique(mysea.mapped)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check unmapped in drugbank for approved/experimental status
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  table(db2rxn$db[`DrugBank ID` %in% mysea.unmapped$intermediate_id]$status.approved)
  ##  approved not_approved 
  ##        9            35 ... ok, so most of the unmapped are likely unapproved drugs, moving on without them ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # save(sc2rxn, file="/sc/arion/projects/psychgen/methods/rx/files/drug_id_map_sc2rxn.Rd")
