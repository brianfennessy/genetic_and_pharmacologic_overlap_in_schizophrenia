# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/full_database.xml
# 2. my.fil = list.files(pattern="DB*")
# 3. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_ALT_IDS.rrf
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/full_database.Rdata
# 2. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/hiddenDBIDmatches.Rd
# 3. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/full_database_IDharmonization.Rdata
# 4. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/dbid_matches.tsv
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  library(data.table)
  library(XML)
  library(foreach)
  library(parallel)
  library(doMC)
  options(cores = parallel:::detectCores())
  registerDoMC(15)

  # Store user home directory and navigate to it
  user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
  setwd(user_home_directory)

  # Get the official home directory of the user
  home_dir <- getwd()

  # Navigate to the raw data zip file and unzip
  setwd(home_dir)
  setwd("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drugbank_id_harmonization")
  
  makemod = function (x) {trimws(gsub(".", " ", gsub("([.])\\1+", "\\1", tolower(make.names(x))), fixed=TRUE), "both")}                           

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dbank full release csv file (no longer an option in June 2021)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # db = fread ("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/all-drug-links.csv", na="")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dbank full release xml file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp_data = xmlParse("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drug_target/drugbank/5.1.1/full_database.xml")
  dbx = xmlToList(tmp_data)
  save(dbx, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drug_target/drugbank/5.1.1/full_database.Rdata")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drug_target/drugbank/5.1.1/full_database.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parse the scraped html files that contain some mappings not in the full drugbank release
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  my.fil = fread("idlist", header=F)
  my.fil = my.fil$V1
  res = foreach(i = 1:length(my.fil), .combine = rbind )%dopar%{
      input = my.fil[i] 
      if (i %% 100 == 0 ) cat("\r",i,"\t\t")
      my.lst = list("dc.identifier"=c(), "dc.title"=c(), "link"=c(), "title"=c())
      tmp = htmlTreeParse(input)
      root = xmlRoot(tmp)
      kids = xmlChildren(root)
      kidA = kids[[1]]
      metx = which(names(xmlChildren(kidA)) == "meta")
      linx = which(names(xmlChildren(kidA)) == "link")
      tilx = which(names(xmlChildren(kidA)) == "title")
      for (j in metx){
          check1 = xmlGetAttr(kidA[[j]], name="name")
          if(!is.null(check1)){
              if (check1=="dc.identifier"){
                  my.lst$dc.identifier = c(my.lst$dc.identifier, xmlGetAttr(kidA[[j]], name="content"))
              } else if (check1=="dc.title") {
                  my.lst$dc.title = c(my.lst$dc.title, xmlGetAttr(kidA[[j]], name="content"))
              }
          }
      }
      for (k in linx){
          check2 = xmlGetAttr(kidA[[k]], name="rel")
          if(!is.null(check2)){
              if (check2=="canonical"){
                  addlink = gsub("https://www.drugbank.ca/drugs/","",xmlGetAttr(kidA[[k]], name="href"))
                  my.lst$link = c(my.lst$link, addlink)
              }
          }
      }
      for (l in tilx){
          my.lst$title = c(my.lst$title,gsub(" - DrugBank","",xmlValue(kidA[[l]])))
      }
      a0 = input
      a1 = paste(my.lst$dc.identifier, collapse="|")
      a2 = paste(my.lst$dc.title, collapse="|")
      a3 = paste(my.lst$link, collapse="|")
      a4 = paste(my.lst$title, collapse="|")
      my.dt = data.table(dbid = a0, dbname=a4, otherid=a1, othername=a2, linkid=a3)
      my.dt
  }
  nrow(res[otherid!=linkid]) #[1] 0 ... so they are redundant
  res[, linkid:=NULL]
  found = res[dbid != otherid & otherid != "",.(dbid,otherid)]
  save(found, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drug_target/drugbank/5.1.1/hiddenDBIDmatches.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drug_target/drugbank/5.1.1/hiddenDBIDmatches.Rd")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# list the unique external databases that drugbank links to in its full release xml file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  slt = c()
  for (i in 1:(length(dbx)-1)){
      cat(i, "\r")
      if (!is.null(dbx[[i]]$`external-identifiers`)){
          newslt = dbx[[i]]$`external-identifiers`
          for(j in 1:length(newslt)){
              slt = c(slt,dbx[[i]]$`external-identifiers`[[j]]$resource)
          }
      }
  }
  slt = unique(slt)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parse the full release xml file so we can map DB IDs to one another
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  identifiers = c("drugbank-id", "cas-number", "unii", "synonyms", "atc-codes", "external-identifiers")
  res = foreach(i = 1:(length(dbx)-1), .combine = rbind )%dopar%{
      if (i %% 100 == 0 ) cat("\r",i,"\t\t")      
      primaryID = dbx[[i]][[1]]$text
      primaryNAME = dbx[[i]]$name
      check.dt = c()
      for (j in identifiers){check.dt = rbind(check.dt,data.table(name = j, indices = which(names(dbx[[i]]) == j)))}
      for(j in 1:nrow(check.dt)){
          index = check.dt[j]$indices
          check.dt[j,isNull := is.null(dbx[[i]][[index]])]
      }
      check.dt.null = check.dt[ isNull==TRUE ]
  
      for (j in unique(check.dt$name)){
  
          ## drugbank ids
          if (j=="drugbank-id"){
              if (j %in% check.dt.null$name) {
                  my.dbids = NA
              } else { 
                  my.idx = check.dt[name==j]$indices
                  if (length(my.idx)==1){
                      my.dbids=NA
                  } else {
                      my.dbids = c()
                      for (id in my.idx[2:length(my.idx)]){my.dbids = c(my.dbids, dbx[[i]][[id]])}
                      my.dbids = paste(my.dbids, collapse="|")
                  }
              }
          }

          ## synonyms            
          if (j=="synonyms"){
              if (j %in% check.dt.null$name) {
  
                  my.syn = NA
              } else { 
                  my.syn = c()
                  for (nl in 1:length(dbx[[i]]$synonyms)){my.syn=c(my.syn,dbx[[i]]$synonyms[[nl]][[1]])}
                my.syn = paste(my.syn, collapse="|")
              }
          }
          
          ## cas
          if (j=="cas-number"){
              if (j %in% check.dt.null$name) {
                  mycas = NA
              } else { 
                  mycas = dbx[[i]]$`cas-number`[[1]]
              }
          }
          
          ## unii
          if (j=="unii"){
              if (j %in% check.dt.null$name) {
                  myuni = NA
              } else { 
                  myuni = dbx[[i]]$`unii`[[1]]
              }
          }
          
          ## atc
          if (j=="atc-codes"){
              if (j %in% check.dt.null$name) {
                  myatc = NA
              } else { 
                  myatc = dbx[[i]]$`atc-codes`$`atc-code`$.attrs
              }
          }
          
          ## external identifiers
          if (j=="external-identifiers"){
              if (j %in% check.dt.null$name) {
                  ei.dt = NA
              } else { 
                  ei.dt = data.table(dbid=primaryID)
                  for (k in 1:length(dbx[[i]]$`external-identifiers`)){
                      col = dbx[[i]]$`external-identifiers`[[k]]$resource
                      val = dbx[[i]]$`external-identifiers`[[k]]$identifier
                      ei.dt[,(col):=val]
                  }
              }
          }
      }

      final = data.table( dbid = primaryID, dbname = primaryNAME, secid = my.dbids, syn = my.syn, cas = mycas, uni = myuni, atc = myatc )
      if (!is.na(ei.dt)) {
          final = merge(final, ei.dt, by="dbid")
      }

      final[, slt[!slt %in% colnames(final)] := NA ]
      final
  }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine the mappings from the scraped html files with mappings from full release xml
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  res = merge(res,found,all.x=TRUE)
  res[ !is.na(otherid), secid := paste(secid,otherid,sep="|")]
  res [, otherid := NULL ] 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how informative are the various IDs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  xxx = c("cas", "uni", "atc", "Drugs Product Database (DPD)", 
          "PubChem Substance", "KEGG Drug", "PharmGKB", "UniProtKB", "Therapeutic Targets Database", 
          "Wikipedia", "ChEMBL", "GenBank", "KEGG Compound", "ChEBI", "PubChem Compound", 
          "ChemSpider", "BindingDB", "IUPHAR", "Guide to Pharmacology", "PDB")
  checker = c()
  for (i in xxx){
      nu = uniqueN(res[ !is.na (get(i)),.(get(i))])
      nt = nrow(res[ !is.na (get(i))])
      fr = nu/nt 
      checker = rbind(checker, data.table( database = i, nUniq = nu, nTot = nt, frac = fr ))
  }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dt = res[, list( syn = unlist(strsplit(syn, "|", fixed=T))), by=list(`dbid`, `dbname`, `secid`, `cas`, `uni`, `atc`, 
                                                                         `Drugs Product Database (DPD)`,`PubChem Substance`, 
                                                                         `KEGG Drug`, `PharmGKB`, `UniProtKB`, `Therapeutic Targets Database`, 
                                                                         `Wikipedia`, `ChEMBL`, `GenBank`, `KEGG Compound`, `ChEBI`, 
                                                                         `PubChem Compound`, `ChemSpider`, `BindingDB`, `IUPHAR`, `Guide to Pharmacology`, `PDB`)]
  dt = dt[, list( secid = unlist(strsplit(secid, "|", fixed=T))), by=list(`dbid`, `dbname`, `cas`, `uni`, `atc`, `syn`, 
                                                                                `Drugs Product Database (DPD)`,`PubChem Substance`,`KEGG Drug`, 
                                                                                `PharmGKB`, `UniProtKB`, `Therapeutic Targets Database`, `Wikipedia`, 
                                                                                `ChEMBL`, `GenBank`, `KEGG Compound`, `ChEBI`, `PubChem Compound`, 
                                                                                `ChemSpider`, `BindingDB`, `IUPHAR`, `Guide to Pharmacology`, `PDB`)]
  dt[, dbname := makemod(dbname)]
  dt[, Wikipedia := makemod(Wikipedia)]
  dt[, syn := makemod(syn)]
  colnames(dt) = tolower(make.names(colnames(dt)))  
  save(dt, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drug_target/drugbank/5.1.1/full_database_IDharmonization.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# map DB IDs to one another
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  idtypes = c("dbname", "cas", "uni", "atc", "syn", "drugs.product.database..dpd.", "pubchem.substance", 
              "kegg.drug", "pharmgkb", "uniprotkb", "therapeutic.targets.database", 
              "wikipedia", "chembl", "genbank", "kegg.compound", "chebi", "pubchem.compound", 
              "chemspider", "bindingdb", "iuphar", "guide.to.pharmacology", "pdb", "secid")

  tracker = foreach(i = 1:nrow(dt), .combine = rbind )%dopar%{
      if (i %% 100 == 0 ) cat("\r",i, "of", nrow(dt),"\t\t")      
      pid = dt[i]$dbid
      newdt = c()
      for (j in idtypes){
          if (!is.na(dt[i,get(j)])){
              oid = dt[i,get(j)]
              newrow = data.table(dbid=pid, othid=oid, othid.type=j)
              newdt = rbind(newdt, newrow)
          }
      }
      newdt
  }
  tracker = unique(tracker[!othid %in% c("NA", "na")])
  tracker[, oth := paste(othid, othid.type)]
  matches = tracker[ nchar(othid)==7 ][grep("DB", othid)][othid.type=="secid",.(id1=dbid, id2=othid, meth="drugbank")]
  matches2 = foreach(i = 1:uniqueN(tracker$dbid), .combine = rbind )%dopar%{
      if (i %% 100 == 0 ) cat("\r",i, "of", uniqueN(tracker$dbid),"\t\t")      
      id = unique(tracker$dbid)[i]
      add = data.table(id1=id, id2=id, meth="self")
      for (j in unique(tracker[dbid==id]$oth)){
          if (length(tracker[oth==j & dbid!=id]$dbid)>=1) {
              add = rbind(add, data.table(id1=id, id2=tracker[oth==j & dbid!=i]$dbid, meth=tracker[oth==j & dbid!=i]$othid.type))
          }
      }
      add
  }
  matches = rbind(matches, matches2)
  matches = matches[!meth %in% c("uniprotkb", "genbank", "self", "drugs.product.database..dpd.", "atc")]
  matches = matches[id1 != id2 ]
  
  given = fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_ALT_IDS.rrf", header=F)[,.(id1=V1, id2=V2, meth="drugbank")] 
  matches = rbind(matches, given)
  write.table( matches, row=F, quo=F, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/drug_target/drugbank/5.1.1/dbid_matches.tsv", sep="\t" ) 
