# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Store user home directory and navigate to it
user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
setwd(user_home_directory)

# Get the official home directory of the user
home_dir <- getwd()

# Navigate to the raw data zip file and unzip
setwd(home_dir)

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
# Load dbid_gene to uniprot_gene mapping
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt1_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_AA_SEQ.rrf"
dt1 <- fread(dt1_path)[,.(db_gene=V9,uniprot=V6)]
dt1 <- unique(dt1)

dt2_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_SYN.rrf"
dt2 <- fread(dt2_path)[,.(db_gene=V5,uniprot=V2)]
dt2 <- unique(dt2)

dt3_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_EXT_IDS.rrf"
dt3 <- fread(dt3_path)[,.(db_gene=V6,uniprot=V3)]
dt3 <- unique(dt3)

dt4_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_EXT_IDS.rrf"
dt4 <- fread(dt4_path)[,.(db_gene=V6,uniprot=V3)]
dt4 <- unique(dt4)

dt5_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_GO_CLASSIFIERS.rrf"
dt5 <- fread(dt5_path)[,.(db_drug=V7,db_gene=V6,uniprot=V3,mechanism=V2)]
dt5 <- unique(dt5)

dt6_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_PFAMS.rrf"
dt6 <- fread(dt6_path)[,.(db_gene=V6,uniprot=V3)]
dt6 <- unique(dt6)

dt7_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGET_POLYPEPTIDES.rrf"
dt7 <- fread(dt7_path, sep="|", fill=TRUE)

dt8_path <- "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_TARGETS_GENE_SEQ.rrf"
dt8 <- fread(dt8_path, sep="|", fill=TRUE)[,.(db_gene=V9,uniprot=V6)]
dt8 <- unique(dt8)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify mapping counts for db_gene and uniprot ID's
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get general pair counts
uniqueN(dt1) # [1] 4985
uniqueN(dt2) # [1] 4540
uniqueN(dt3) # [1] 4985
uniqueN(dt4) # [1] 4985
uniqueN(dt5[,.(db_gene,uniprot)]) # [1] 4941
uniqueN(dt6) # [1] 4914
uniqueN(dt7) # Not a good file to use
uniqueN(dt8) # [1] 5255

# See if we have any NA or space values and get rid of those we do
dt1[db_gene=="" | uniprot=="" | is.na(db_gene) | is.na(uniprot),] # Empty data.table (0 rows and 2 cols): db_gene,uniprot
dt2[db_gene=="" | uniprot=="" | is.na(db_gene) | is.na(uniprot),] # Empty data.table (0 rows and 2 cols): db_gene,uniprot
dt3[db_gene=="" | uniprot=="" | is.na(db_gene) | is.na(uniprot),] # Empty data.table (0 rows and 2 cols): db_gene,uniprot
dt4[db_gene=="" | uniprot=="" | is.na(db_gene) | is.na(uniprot),] # Empty data.table (0 rows and 2 cols): db_gene,uniprot
dt5[db_gene=="" | uniprot=="" | is.na(db_gene) | is.na(uniprot),] # Empty data.table (0 rows and 2 cols): db_gene,uniprot
dt6[db_gene=="" | uniprot=="" | is.na(db_gene) | is.na(uniprot),] # Empty data.table (0 rows and 2 cols): db_gene,uniprot
dt8[db_gene=="" | uniprot=="" | is.na(db_gene) | is.na(uniprot),]
#     db_gene uniprot
#   1:             247
#   2:            1309
#   3:            1351
#   4:            1455
#   5:            1507
#  ---                
# 499:           18202
# 500:           18203
# 501:           18427
# 502:           18446
# 503:    
dt8 <- dt8[db_gene!="" & uniprot!="" & !is.na(db_gene) & !is.na(uniprot),]
uniqueN(dt8) # [1] 4752

# Isolate dt5 so we only have the mapping information for now
# Keep in mind that this could be a nice file to refer back to (as well as 7) for more lenghty descriptions!
dt5 <- unique(dt5[,.(db_gene,uniprot)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rbind all files to create a master list of pairs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final <- unique(rbind(dt1,dt2,dt3,dt4,dt5,dt6,dt8))

uniqueN(final$db_gene) # [1] 4761
uniqueN(final$uniprot) # [1] 4766
uniqueN(final[,.(db_gene, uniprot)]) # [1] 4985

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See what kind of an overlap we have between the gene dbid's from our mapping files and our actions file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
act <- unique(actions$dbid_gene)
dt <- unique(final$db_gene)

length(act) # [1] 4865
length(dt) # [1] 4761
length(intersect(act,dt)) # [1] 4761
4761/4865 # [1] 0.9786228 --> Looks like we have uniprot codes for 4761 out of 4865 DB gene ID's, or 97.9%, of DB gene ID's

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in drugbank mapping file for uniprot to gene name
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dbkmap <- fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drugbank_uniprot_id_map.tsv", na="", col.names = c("uniprot","gene","gene_other")) 

uni <- unique(dbkmap$uniprot)
db <- unique(final$uniprot)

length(uni) # [1] 4981
length(db) # [1] 4766
length(intersect(uni,db)) # [1] 4757
4757/4766 # [1] 0.9981116 --> Looks like we have gene names for 4757 out of 4766 uniprot codes , or 99.8%, of uniprot codes

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge by uniprot code to get gene names in out actions file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final <- merge(final,dbkmap,by="uniprot",all.x=TRUE)
uniqueN(final[is.na(gene),]) # [1] 382
uniqueN(final) # [1] 4985
382/4985 # [1] 0.07662989 --> We lose 7.7% of our set when adding gene name

final <- final[!is.na(gene)]
final[,gene:=toupper(gene)]
final[,gene_other:=toupper(gene_other)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(final, sep="\t", "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_target/drugbank/5.1.1/dbgeneid2uni.tsv")
