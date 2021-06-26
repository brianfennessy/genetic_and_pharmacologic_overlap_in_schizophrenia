# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ ME
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/README
###---------------------------------------------------------------------------------
### NOTE on stitch ids used in sider (from http://stitch3.embl.de/download/README): 
###---------------------------------------------------------------------------------
### 
### Chemicals are derived from PubChem. As described in the STITCH paper,
### we merge salt forms and isomers. However, since STITCH 3, isomeric
### compounds can also be investigated separately. In the download files,
### the following convention holds:
###
### CID0... - this is a stereo-specific compound, and the suffix is the PubChem compound id.
###
### CID1... - this is a "flat" compound, i.e. with merged stereo-isomers, and the suffix (without the leading "1") is the PubChem compound id.
###
###----------------------------------------------------
### AWC addendum to above note: 
###----------------------------------------------------
### CID0 is used in stitch version 1 and 2 for the flat IDs
###
###----------------------------------------------------
### SIDER FILE COLUMNS
###----------------------------------------------------
### for full details /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/README
###
###----------------------------------------------------
### NOTE ON UMLS IDs IN SIDER (IN meddra_all_se.tsv)
###----------------------------------------------------
###
###  Sider gives 2 umls ids
###
###       1. "UMLS concept id as it was found on the label"
###       2. "UMLS concept id for MedDRA term"    <--------------- I USED THIS ONE
###                            
###       See email from Michael Kuhn (sider developer) on 1/4/17 that explains the difference.                      
###                            
###  Sider gives "MedDRA concept type" as well, but I excluded because I found that it made things 
###       uneccessarily confusing. Instead, I am just using the UMLS meddra id and the meddra
###       concept name that sider gives and mapping it to the official meddra release myself. 
###
###----------------------------------------------------
### NOTE ON AE FREQ IN SIDER
###----------------------------------------------------
###  In many instances, multiple frequency ranges are presented for same drug-ae (due to multiple labels)
###  I removed the "description" column from freq file, and just took the mean lower and upper bounds
###  so that for any drug-ae pair there is only 1 row.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/stitch/chemical.sources.v4.0.tsv
# 2. /sc/arion/projects/psychgen/methods/rx/data/stitch/chemicals.inchikeys.v4.0.1.tsv
# 3. /sc/arion/projects/psychgen/methods/rx/data/stitch/pc2atc.txt
# 4. /sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata
# 5. /sc/arion/projects/psychgen/methods/rx/data/uniII/UNIIs_10Nov2016_Records.txt
# 6. /sc/arion/projects/psychgen/methods/rx/data/umls/2016AB/installation1/2016AB/META/MRFILES.RRF
# 7. /sc/arion/projects/psychgen/methods/rx/data/umls/2016AB/installation1/2016AB/META/MRCONSO.RRF
# 8. /sc/arion/projects/psychgen/methods/rx/data/umls/2016AB/installation1/2016AB/META/MRSAB.RRF
# 9. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/sider_awcQC.txt
# 10. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/sider_ind_awcQC.txt
# 11. /sc/arion/projects/psychgen/methods/rx/data/stitch/src14src22.txt
# 12. /sc/arion/projects/psychgen/methods/rx/data/stitch/pc2unii.txt
# 13. /sc/arion/projects/psychgen/methods/rx/data/pcie/input/siderPC_st4.txt
# 14. /sc/arion/projects/psychgen/methods/rx/data/pcie/output/OUT_siderPC.txt
# 15. /sc/arion/projects/psychgen/methods/rx/files/pc_compound_type.ttl.reformatted
# 16. /sc/arion/projects/psychgen/methods/rx/data/pcie/input/rxn2smil2pc.txt
# 17. /sc/arion/projects/psychgen/methods/rx/data/pcie/output/OUT_rxn2smil2pc.txt
# 18. /sc/arion/projects/psychgen/methods/rx/data/stitch/chemicals.v4.0.tsv
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_sdr2rxn.Rd
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
  setwd("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/")
  scratchdir = "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/"
  library(interval)
  library(bit)
  library(bit64)
  library(stringr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# stitch 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stitch4.fil = './genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/stitch/chemical.sources.v4.0.tsv'
  stitch4.fl2 = './genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/stitch/chemicals.inchikeys.v4.0.1.tsv'
  stitch4.src = fread(stitch4.fil, col.names=c("stitchF", "stitchS","source","source_id"),fill=TRUE, na=c("","NA")) 
  stitch4.chi = fread(stitch4.fl2, header=T, sep='\t') 
  pc2sti = stitch4.src [ source=="PC",  c(1,2,4),with=F ]
  colnames(stitch4.chi)[c(1,2,3)] = colnames(pc2sti)
  unique(stitch4.src$source) #shows databases stitch links to
  colnames(pc2sti)[3]="pubchem"
  pc2sti[,pubchem := as.numeric(pubchem)]
  pc2atcST = stitch4.src [ source=="ATC", c(1,2,4),with=F ]
  colnames(pc2atcST)[3]="atc"
  pc2atcST = merge(pc2sti, pc2atcST, by=c("stitchF", "stitchS"))[,c("atc","pubchem"),with=F]
  pc2atcPC = fread("../../stitch/pc2atc.txt", col.names=c("pubchem","atc"))[,c("atc","pubchem"),with=F]
  pc2atc = unique(rbind(pc2atcPC,pc2atcST))
  pc2atc2sti = merge(pc2sti, pc2atc, by="pubchem")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rxnorm 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata")  
  rxnLatest.conso[, STR.mod:=makemod(STR)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# unii 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  unii = fread('./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/uniII/UNIIs_10Nov2016_Records.txt', sep='\t', na.strings=c("NA","") )
  unii[, PT.mod := makemod(PT)]
  unii[,RN.mod := str_pad(gsub("-","",RN), 10, pad="0")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# umls 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uml = "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/umls/2016AB/installation1/2016AB/META/"
  uml.cols = fread(paste0(uml,"MRFILES.RRF"), header=F, sep="|")
  uml.con = fread(paste0(uml,'MRCONSO.RRF'), sep="|",header=F)
  uml.sab = fread(paste0(uml,'MRSAB.RRF'), sep="|",header=F)[,.(sab=V4,sab_full=V5)]
  uml.con[,V19:=NULL]
  colnames(uml.con) = unlist(strsplit(uml.cols[V1=="MRCONSO.RRF"]$V3, split=","))
  uml.con[, STR.mod := makemod(STR)]
  rxn2uml = unique(uml.con[ SAB=="RXNORM",c("CODE","CUI"),with=F])
  colnames(rxn2uml)=c("rxcui","umlcui")
  ## merge(as.data.table(table(uml.con[CUI %in% uml.con[SAB=="RXNORM",]$CUI & SAB!="RXNORM",]$SAB)),uml.sab,by.x="V1",by.y="sab") #[ N > 10,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdr.adv = fread("sider_awcQC.txt", sep='\t', header=T)
  sdr.ind = fread("sider_ind_awcQC.txt", header=T, sep='\t')
  sdr.ind = unique(merge(sdr.ind, stitch4.src, by.x="sider.stitchF", by.y="stitchF")[,.(sider.stitchF,sider.stitchS=stitchS,sider.meddra_name_umls,meddra,lev)])
  sdr.pairs = unique(rbind(sdr.adv[,.(sider.stitchF,sider.stitchS)], sdr.ind[,.(sider.stitchF,sider.stitchS)]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# stitch4 subset for sider drugs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stitch4.src.sdr = merge( stitch4.src, sdr.pairs, by.x = c("stitchF", "stitchS"), by.y=c("sider.stitchF", "sider.stitchS"))
  stitch4.src.sdr.mapped = c()
  stitch4.src.sdr.unmapped = copy(stitch4.src.sdr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider <> rxnorm, via drugbank
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stitch4.src.sdr.db = unique(stitch4.src.sdr[ source=="DrugBank",.(stitchF, dbank=source_id) ])
  rxn2db = unique(rxnLatest.conso[SAB=="DRUGBANK",.(RXCUI, dbank=CODE)])
  stitch2rxn.db = unique(merge( stitch4.src.sdr.db, rxn2db )[,.(stitchF,RXCUI)])
  found = get_rxcui_in(stitch2rxn.db$RXCUI, scratchdir, use.foreach=FALSE)
  stitch2rxn.db = unique(merge( stitch2rxn.db[,.(sider.stitchF=stitchF, rxcui_input=RXCUI)], found, by="rxcui_input")[,.(sider.stitchF,rxcui_output)])
  stitch4.src.sdr.mapped = rbind(stitch4.src.sdr.mapped, stitch2rxn.db)
  stitch4.src.sdr.unmapped = stitch4.src.sdr.unmapped[ ! stitchF %in% stitch4.src.sdr.mapped$sider.stitchF ]
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) / uniqueN(stitch4.src.sdr$stitchF) #[1] 0.6443265 ... verified 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider <> rxnorm, via atc
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stitch4.src.sdr.atc = unique(stitch4.src.sdr.unmapped[ source=="ATC",.(stitchF, atc=source_id) ])
  rxn2atc = unique(rxnLatest.conso[ SAB=="ATC" & CODE != "NOCODE", .(RXCUI, atc=CODE)] )
  stitch2rxn.atc = unique(merge( stitch4.src.sdr.atc, rxn2atc )[,.(stitchF,RXCUI)])
  found = get_rxcui_in(stitch2rxn.atc$RXCUI, scratchdir, use.foreach=FALSE)
  stitch2rxn.atc = unique(merge( stitch2rxn.atc[,.(sider.stitchF=stitchF, rxcui_input=RXCUI)], found, by="rxcui_input")[,.(sider.stitchF,rxcui_output)])
  stitch4.src.sdr.mapped = rbind(stitch4.src.sdr.mapped, stitch2rxn.atc)
  stitch4.src.sdr.unmapped = stitch4.src.sdr.unmapped[ ! stitchF %in% stitch4.src.sdr.mapped$sider.stitchF ]
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) / uniqueN(stitch4.src.sdr$stitchF) #[1] 0.7803583 ... in 13apr2021 is 0.7810219 ...

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider <> rxnorm, via unii
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rxn2unii = unique( unii[which(!is.na(unii$RXCUI))   ,c("UNII","RXCUI"),with=F]);colnames(rxn2unii) = c("unii","rxcui") # rxn <> unii
  pc2uniiUC = fread('./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/stitch/src14src22.txt', header=T, col.names = c("unii","pubchem"), sep='\t' )
  pc2uniiPC = fread('./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/stitch/pc2unii.txt', header=T, col.names = c("unii","pubchem"), sep='\t' )
  pc2unii = unique(rbind(pc2uniiUC, pc2uniiPC))
  pc2unii2sti = merge(pc2sti, pc2unii, by="pubchem")[stitchF %in% stitch4.src.sdr.unmapped$stitchF]
  stitch2rxn.uni = merge(pc2unii2sti, rxn2unii, by="unii")
  found = get_rxcui_in(stitch2rxn.uni$rxcui, scratchdir, use.foreach=FALSE)
  stitch2rxn.uni = unique(merge( stitch2rxn.uni[,.(sider.stitchF=stitchF, rxcui_input=rxcui)], found, by="rxcui_input")[,.(sider.stitchF,rxcui_output)])
  stitch4.src.sdr.mapped = rbind(stitch4.src.sdr.mapped, stitch2rxn.uni)
  stitch4.src.sdr.unmapped = stitch4.src.sdr.unmapped[ ! stitchF %in% stitch4.src.sdr.mapped$sider.stitchF ]
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) / uniqueN(stitch4.src.sdr$stitchF) #[1] 0.9216987 ... in 13apr2021 is 0.9223623

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider <> rxnorm, via spl
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## rxnorm <> spl
  uml2spl = unique(uml.con[CUI %in% unique(rxn2uml$umlcui) & SAB=="MTHSPL",c("CUI","CODE"),with=F])
  colnames(uml2spl)=c("umlcui","mthspl")
  rxn2splUM = unique(merge( rxn2uml, uml2spl, by="umlcui" )[,c("rxcui","mthspl"),with=F])
  rxn2splRX = unique(rxnLatest.conso[ SAB=="MTHSPL" & CODE != "NOCODE", c("RXCUI","CODE"),with=F])
  colnames(rxn2splRX) = c("rxcui","mthspl")
  rxn2spl = unique(rbind(rxn2splUM,rxn2splRX))
  
  ## pubchem <> spl
  ##make pcie input file (settings on pcie: input id list = CIDs; output IDs = fda/spl indexing data; operator type = same cid)
  pcie_in = './genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/pcie/input/siderPC_st4.txt'
  pcie_out = './genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/pcie/output/OUT_siderPC.txt'
  write.table(unique(stitch4.src.sdr[ source=="PC","source_id", ]), row.names=F, col.names=F, quote=F, file = pcie_in) #(read in pcie output)
  tmp = unique(fread(pcie_out, sep='\t', header=F, col.names=c("pubchem","mthspl") ) [mthspl != "",])
  pc2spl2sti  = merge(pc2sti, tmp, by="pubchem" )[stitchF %in% stitch4.src.sdr.unmapped$stitchF]
  
  ## add
  stitch2rxn.spl = merge(pc2spl2sti, rxn2spl, by="mthspl")
  found = get_rxcui_in(stitch2rxn.spl$rxcui, scratchdir,use.foreach=F) #use.foreach=F if on node without internet
  stitch2rxn.spl = unique(merge( stitch2rxn.spl[,.(sider.stitchF=stitchF, rxcui_input=rxcui)], found, by="rxcui_input")[,.(sider.stitchF,rxcui_output)])
  stitch4.src.sdr.mapped = rbind(stitch4.src.sdr.mapped, stitch2rxn.spl)
  stitch4.src.sdr.unmapped = stitch4.src.sdr.unmapped[ ! stitchF %in% stitch4.src.sdr.mapped$sider.stitchF ]
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) / uniqueN(stitch4.src.sdr$stitchF) #[1] 0.9270073

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider <> rxnorm, via ndf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## rxnorm <> ndf
  uml2ndf = unique(uml.con[CUI %in% unique(rxn2uml$umlcui) & SAB=="NDFRT" ,c("CUI","CODE")])
  colnames(uml2ndf)=c("umlcui","ndfrt")
  rxn2ndfUM = unique(merge( rxn2uml, uml2ndf, by="umlcui" )[,c("rxcui","ndfrt"),with=F])
  rxn2ndfRX = unique(rxnLatest.conso[ SAB=="NDFRT"  & CODE != "NOCODE", c("RXCUI","CODE"),with=F])
  colnames(rxn2ndfRX) = c("rxcui","ndfrt")
  rxn2ndf = unique(rbind(rxn2ndfUM,rxn2ndfRX))

  ## pubchem <> ndf
  pc2ndf = fread('./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/pc_compound_type.ttl.reformatted', sep='\t', header=F) [V1!="",c(1,3),with=F]
  pc2ndf$V3 = do.call('rbind', strsplit(as.character(pc2ndf$V3), ' ', fixed=T))[,1]
  pc2ndf$V1 = do.call('rbind', strsplit(as.character(pc2ndf$V1), 'CID', fixed=T))[,2]
  pc2ndf = pc2ndf[grep("ns4",pc2ndf$V3),]
  pc2ndf$V3 = do.call('rbind', strsplit(as.character(pc2ndf$V3), ':', fixed=T))[,2]
  colnames(pc2ndf) = c("pubchem","ndfrt")
  pc2ndf[,pubchem:=as.numeric(pubchem)]
  pc2ndf2sti = merge(pc2sti, pc2ndf, by="pubchem")[stitchF %in% stitch4.src.sdr.unmapped$stitchF]

  ## final map
  stitch2rxn.ndf = merge(pc2ndf2sti, rxn2ndf, by="ndfrt")
  found = get_rxcui_in(stitch2rxn.ndf$rxcui, scratchdir, use.foreach=F) 
  stitch2rxn.ndf = unique(merge( stitch2rxn.ndf[,.(sider.stitchF=stitchF, rxcui_input=rxcui)], found, 
                                by="rxcui_input", allow.cartesian=TRUE) [,.(sider.stitchF,rxcui_output)])
  stitch4.src.sdr.mapped = rbind(stitch4.src.sdr.mapped, stitch2rxn.ndf)
  stitch4.src.sdr.unmapped = stitch4.src.sdr.unmapped[ ! stitchF %in% stitch4.src.sdr.mapped$sider.stitchF ]
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) / uniqueN(stitch4.src.sdr$stitchF) #[1] 0.9329794 ... 13apr2021 is 0.9296616

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider <> rxnorm, via smiles
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## rxn <> smiles
  pcieR_in = './genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/pcie/input/rxn2smil2pc.txt'
  pcieR_out = './genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/pcie/output/OUT_rxn2smil2pc.txt'
  rxn2smi  = unique( unii[which(!is.na(unii$RXCUI) & !is.na(unii$SMILES))  ,c("RXCUI","SMILES"),with=F])
  write.table(unique(rxn2smi[,"SMILES",with=F]), row.names=F, col.names=F, quote=F, file=pcieR_in)
  tmp = unique(fread(pcieR_out, sep='\t', header=F, col.names=c("smiles","pubchem") ) [!is.na(pubchem),])

  ## pubchem <> smiles
  pc2rxn.bysmi = unique(merge( tmp, rxn2smi, by.y="SMILES",by.x="smiles")[,.(pubchem,rxcui=RXCUI)])
  stitch2rxn.smi = merge(pc2rxn.bysmi, pc2sti, by="pubchem")[stitchF %in% stitch4.src.sdr.unmapped$stitchF]
  found = get_rxcui_in(stitch2rxn.smi$rxcui, scratchdir,use.foreach=F) #use.foreach=F if on node without internet
  stitch2rxn.smi = unique(merge( stitch2rxn.smi[,.(sider.stitchF=stitchF, rxcui_input=rxcui)], found, 
                                by="rxcui_input", allow.cartesian=TRUE) [,.(sider.stitchF,rxcui_output)])
  stitch4.src.sdr.mapped = rbind(stitch4.src.sdr.mapped, stitch2rxn.smi)
  stitch4.src.sdr.unmapped = stitch4.src.sdr.unmapped[ ! stitchF %in% stitch4.src.sdr.mapped$sider.stitchF ]
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) / uniqueN(stitch4.src.sdr$stitchF) #0.9402787 ... 13apr2021 its 0.938288

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider <> rxnorm, via inchikey
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rxn2chi  = unique( unii[!is.na(RXCUI) & !is.na(INCHIKEY),.(rxcui=RXCUI,inchikey=INCHIKEY)]) ## rxnorm <> inchkey
  pc2chi2sti = stitch4.chi # pubchem <> inchikey
  ##pc2chi = unique(pc2chi2sti[,.(pubchem, inchikey)])
  pc2rxn2sti.bychi = merge( pc2chi2sti, rxn2chi, by="inchikey")[stitchF %in% stitch4.src.sdr.unmapped$stitchF]

  ## how much will this add?
  uniqueN(pc2rxn2sti.bychi$stitchF)/uniqueN(stitch4.src.sdr$stitchF) ##[1] 0 ... still zero in 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# remaining unmapped sider terms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stitch4 = fread('./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/stitch/chemicals.v4.0.tsv', header=T, sep='\t')
  unmapped.st4 = unique(stitch4[chemical %in% c(stitch4.src.sdr.unmapped$stitchF,stitch4.src.sdr.unmapped$stitchS),.(chemical, name=makemod(name))])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider <> rxnormn, by string
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stitch2rxn.str = merge(unmapped.st4, unique(rxnLatest.conso[,.(RXCUI, name=STR.mod)])) 
  found = get_rxcui_in(unique(stitch2rxn.str$RXCUI), scratchdir,use.foreach=F) 
  stitch2rxn.str = unique(merge( stitch2rxn.str[,.(sider.stitchF=chemical, rxcui_input=RXCUI)], found, 
                                by="rxcui_input", allow.cartesian=TRUE) [,.(sider.stitchF,rxcui_output)])
  stitch2rxn.strF = stitch2rxn.str[grep("CID1", sider.stitchF)]
  stitch2rxn.strS = stitch2rxn.str[grep("CID0", sider.stitchF)]
  stitch2rxn.strS = unique(merge(stitch2rxn.strS, stitch4.src.sdr.unmapped, by.x="sider.stitchF", by.y="stitchS")[,.(sider.stitchF=stitchF, rxcui_output)])
  stitch2rxn.str = unique(rbind(stitch2rxn.strF, stitch2rxn.strS))  
  stitch4.src.sdr.mapped = rbind(stitch4.src.sdr.mapped, stitch2rxn.str)
  stitch4.src.sdr.unmapped = stitch4.src.sdr.unmapped[ ! stitchF %in% stitch4.src.sdr.mapped$sider.stitchF ]
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) / uniqueN(stitch4.src.sdr$stitchF) #0.9615129 ... confirmed same in 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# still unmapped 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  unmapped.st4 = unique(stitch4[chemical %in% c(stitch4.src.sdr.unmapped$stitchF,stitch4.src.sdr.unmapped$stitchS),.(chemical, name=makemod(name))])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check number of the 1507 flat stitch IDs in sider have mapped to >= 1 RXCUI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) #[1] 1449 ... confirmed same in 13apr2021
  uniqueN(stitch4.src.sdr.mapped$sider.stitchF) / uniqueN(stitch4.src.sdr$stitchF) #[1] 0.9615129 ... confirmed same in 13apr2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdr2rxn = copy(stitch4.src.sdr.mapped)
  save(sdr2rxn, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_sdr2rxn.Rd")
