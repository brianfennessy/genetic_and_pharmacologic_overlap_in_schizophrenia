# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ ME
# 1. Drugbank has its own ID system (eg "DB12345")
# 2. RXNORM has integrated Drugbank as a vocabulary (ie, SAB as "DRUGBANK")
# 3. Many of the DB ids in RXNORM however do not map to RXCUI that has a SAB value of "RXNORM". 
#    Only such RXCUIs are in the relationship file. 
#    Therefore, a RXCUI that doesn't map to an RXCUI with SAB of RXNORM is basically an unmapped term.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata
# 2. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/all-drug-links.csv
# 3. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_GROUPS.rrf
# 4. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/hiddenDBIDmatches.Rd
# 5. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_ALT_IDS.rrf
# 6. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_SYN.rrf
# 7. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_PROD.rrf
# 8. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_INT_BRANDS.rrf
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/db2rxn.tsv
# 2. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_db2rxn.Rd
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  makemod = function (x) {trimws(gsub(".", " ", gsub("([.])\\1+", "\\1", tolower(make.names(x))), fixed=TRUE), "both")}                           
  library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in latest rxnorm (9/4/18 release)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("/sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata")
  rxnLatest.conso[, STR.mod := makemod(STR) ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in dbank
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  db = fread ("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/all-drug-links.csv", na="")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add approval status to dbank data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  a0 = fread("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_GROUPS.rrf", sep="|",na="", header=F)
  a0 = a0[, list(status=paste(V1, collapse="|")), by=list(V2) ]
  db = merge(db, a0, by.x="DrugBank ID", by.y="V2") 
  db[, rxnmapped.any := "not.mapped.to.any.sab" ]
  db[, rxnmapped.rxn := "not.mapped.to.rxn.sab" ]
  db[, status.approved := "not_approved" ]
  db[grep("approved", status), status.approved := "approved" ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of dbids
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dbids = unique(db$`DrugBank ID`)
  length( dbids ) #[1] 11292 ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of dbids in rxnorm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dbids.inrxn = unique( rxnLatest.conso[CODE %in% dbids]$CODE )
  length( dbids.inrxn ) #6335 ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of rxcui linked to a dbid
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  db[ `DrugBank ID` %in% dbids.inrxn, rxnmapped.any:="mapped.to.any.sab"]
  dbids.rxcui = unique( rxnLatest.conso[CODE %in% dbids]$RXCUI )
  length(dbids.rxcui) #[1] 8007 ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of dbids linked to rxcui with sab == "rxnorm"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dbids.oksab = unique( rxnLatest.conso[RXCUI %in% dbids.rxcui & SAB=="RXNORM"]$RXCUI )
  dbids.oksab = unique( rxnLatest.conso[RXCUI %in% dbids.oksab & SAB=="DRUGBANK" & CODE %in% dbids]$CODE )
  length(dbids.oksab) #[1] 3617 ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of dbids linked to rxcui with sab != "rxnorm"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  db[ `DrugBank ID` %in% dbids.oksab, rxnmapped.rxn:="mapped.to.rxn.sab"]
  dbids.dbsab = dbids.inrxn[!dbids.inrxn %in% dbids.oksab]
  length(dbids.dbsab) #[1] 2718 ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of dbids either not in rxnorm OR in rxnorm with sab != "rxnorm"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dbids.unmapped = unique(c(dbids.dbsab, dbids[!dbids %in% dbids.inrxn]))
  length(dbids.unmapped) #[1] 7675 ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sanity checks
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  2718 + 3617 #[1] 6335 ... ok good ... confirmed 13apr2021
  7675 + 3617 #[1] 11292 ...ok good ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge dbank and rxn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  db.mapped = db[rxnmapped.rxn=="mapped.to.rxn.sab",.(dbid=`DrugBank ID`)]
  db.mapped = merge(db.mapped, rxnLatest.conso[,.(dbid=CODE, rxcui=RXCUI)], by="dbid")
  db.mapped = unique(merge(db.mapped, rxnLatest.conso[SAB=="RXNORM",.(rxcui=RXCUI)], by="rxcui"))
  # db.mapped = db.mapped[,list(rxcui=paste(rxcui,collapse="|")), by=list(dbid)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in data mapping DBIDs to one another
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/hiddenDBIDmatches.Rd")
  a0 = fread("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_ALT_IDS.rrf", sep="|",na="", header=F)
  found = rbind(found, a0[nchar(V1)==7][grep("^DB", V1)][grep("^DB",V2),.(dbid=V1, otherid=V2)])
  db[, `Mapped DrugBank ID`:=`DrugBank ID`]
  for (i in 1:nrow(found)) { 
      id1 = found[i]$dbid
      id2 = found[i]$otherid
      status1 = db[`DrugBank ID` == id1]$rxnmapped.rxn
      status2 = db[`DrugBank ID` == id2]$rxnmapped.rxn
      if (length(status1)>0 & length(status2)>0){
          if (status1=="mapped.to.rxn.sab" & status2=="not.mapped.to.rxn.sab"){
              cat(i, "CHANGING\n")
              db[`DrugBank ID` == id2, rxnmapped.rxn:="mapped.to.rxn.sab"]
              db[`DrugBank ID` == id2, `Mapped DrugBank ID`:=id1]
          } else if (status2=="mapped.to.rxn.sab" & status1=="not.mapped.to.rxn.sab"){
              cat(i, "CHANGING\n")
              db[`DrugBank ID` == id1, rxnmapped.rxn:="mapped.to.rxn.sab"]
              db[`DrugBank ID` == id1, `Mapped DrugBank ID`:=id2]
          } else {
              cat(i,"\n")
          }
      } else {
          cat(i,"\n")
      }
  }
  new = db[`DrugBank ID` != `Mapped DrugBank ID`,.(dbid1=`DrugBank ID`, dbid2=`Mapped DrugBank ID`)]
  db.mapped = rbind(db.mapped, merge(new, db.mapped, by.x="dbid2", by.y="dbid")[,.(dbid=dbid1, rxcui)]) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# update unmapped list
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dbids.unmapped = db[!`DrugBank ID` %in% db.mapped$dbid]$`DrugBank ID`

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# secondary strings matching to rxnorm 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  a1 = fread("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_SYN.rrf", sep="|",na="", header=F)
  a1[, V1:= gsub("'", "", V1)]
  a1 = a1[, .(dbid=V4, term=makemod(V1), termtype="syn")]
  a2 = fread("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_PROD.rrf", sep="|",na="", header=F)
  a2 = a2[, .(dbid=V16, term=makemod(V1), termtype="prod")]
  a3 = fread("/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/DRUGBANK_DRUG_INT_BRANDS.rrf", sep="|",na="", header=F)
  a3 = a3[, .(dbid=V3, term=makemod(V1), termtype="intbrand")]
  a4 = db[,.(dbid=`DrugBank ID`, term=makemod(Name), termtype="name") ]
  a5 = db[!is.na(`Wikipedia ID`),.(dbid=`DrugBank ID`, term=makemod(`Wikipedia ID`), termtype="wiki") ]
  aa = rbind( a1, a2, a3, a4, a5 )[dbid %in% dbids.unmapped ]
  xxx = merge(aa, rxnLatest.conso[,.(term=STR.mod, rxcui=RXCUI)], by="term")
  xxx = merge(xxx[,.(dbid,rxcui)], rxnLatest.conso[SAB=="RXNORM",.(term=STR.mod, rxcui=RXCUI)], by=c("rxcui"))
  # xxx = unique(xxx)[,list(rxcui=paste(rxcui,collapse="|")), by=list(dbid)]
  xxx = unique(xxx[,.(rxcui,dbid)])
  db.mapped = rbind(db.mapped, xxx)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how did we do
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  nrow(db[`DrugBank ID` %in% db.mapped$dbid & status.approved == "approved"]) / nrow(db[status.approved=="approved"]) #[1] 0.9336678 ... so we map 93% of approved drugs ... confirmed 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # write.table(db.mapped, quo=F, row=F, sep="\t", file="/sc/arion/projects/psychgen/methods/rx/data/drug_target/drugbank/5.1.1/db2rxn.tsv")
  my.db = list("db"=db, "db.mapped"=db.mapped) 
  db2rxn = copy(my.db)
  # save(db2rxn, file="/sc/arion/projects/psychgen/methods/rx/files/drug_id_map_db2rxn.Rd")
