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
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/meddra_all_se.tsv
# 2. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/meddra_freq.tsv
# 3. mdr = parsemdr(version="v20.0")
# 4. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/meddra.tsv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/sider_awcQC.txt
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SIDER v4.1 uses meddra 16.1 
### (per http://sideeffects.embl.de/ - "The current version (SIDER 4.1) has been released on October 21, 2015. This release uses the MedDRA dictionary (version 16.1)")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  setwd("/sc/arion/projects/psychgen/methods/rx/data/side_effects")
  source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in sider data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sidA = unique(fread('sider4.1/meddra_all_se.tsv', header=F, sep='\t') [,c(1,2,5,6),with=F])
  sidB = unique(fread('sider4.1/meddra_freq.tsv',   header=F, sep='\t') [,c(1,2,6,7,9,10),with=F])
  colnames(sidA) = c("stitchF","stitchS","umls_meddra","ae_name")
  colnames(sidB) = c("stitchF","stitchS","freq_lowBound","freq_upBound","umls_meddra","ae_name")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calc mean freq vals for drug/ae pairs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmplo = aggregate(sidB[,"freq_lowBound",drop=F], by=list(sidB$stitchF,sidB$stitchS,sidB$umls_meddra,sidB$ae_name),mean)
  tmphi = aggregate(sidB[,"freq_upBound",drop=F], by=list(sidB$stitchF,sidB$stitchS,sidB$umls_meddra,sidB$ae_name),mean)
  colnames(tmplo) = c("stitchF", "stitchS", "umls_meddra", "ae_name","freq_lowBound")
  colnames(tmphi) = c("stitchF", "stitchS", "umls_meddra", "ae_name", "freq_upBound")
  sidB = as.data.table( merge ( tmplo,tmphi,by=c("stitchF", "stitchS", "umls_meddra", "ae_name") ) )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine ae freq tables into 1 (unclear why released as separate tables)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdr = merge (sidA,sidB,by=colnames(sidA),all.x=T) #there are 3 rows in freq file not in main file, excluding (ie, drug-ae pairs where a frequcny is reported but its not in the main file)
  colnames(sdr) = c("sider.stitchF", "sider.stitchS", "sider.meddra_name_umls", "sider.meddra_name", "sider.freqlo","sider.freqhi")
  sdr = sdr[sider.meddra_name_umls!=""]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in meddra data from official meddra release
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdr = parsemdr(version="v20.0")
  mdr.wide = copy(mdr)
  mdr = rbind(mdr[,.(meddra=llt_code, name=tolower(llt_name),  meddra_level="llt")],
              mdr[,.(meddra=pt_code,  name=tolower(pt_name),   meddra_level="pt")],
              mdr[,.(meddra=hlt_code, name=tolower(hlt_name),  meddra_level="hlt")],
              mdr[,.(meddra=hlgt_code,name=tolower(hlgt_name), meddra_level="hlgt")],
              mdr[,.(meddra=soc_code, name=tolower(soc_name),  meddra_level="soc")])
  mdr[,name_parse := makemod(name)]
  mdr = unique(mdr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in meddra data from sider 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdrmdr = unique(fread('sider4.1/meddra.tsv', header=F, sep='\t')[,.(sider.meddra_name_umls=V1,meddra=V3)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check fraction of sider meddra terms (v16) that map to meddra v20
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sum(sdrmdr$meddra %in% mdr$meddra) / length(sdrmdr$meddra) #[1] 0.9967994
  length(sdrmdr$meddra) - sum(sdrmdr$meddra %in% mdr$meddra) #[1] 238
  ## ... ok, so we only lose ~200 ... (to map between versions to keep there need sinai pw and login for the MVAT tool, which changes each year (contact levy))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# map meddra terms supplied by sider (meddra v16) to meddra v20 pt level (and above)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdrmdr = merge(mdr, sdrmdr)
  lev = sdrmdr[, list( levels = paste(unique(meddra_level), collapse="|")), by=list(sider.meddra_name_umls) ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for terms that are supplied as llt level only, map to pt
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  lev.llt = lev[levels=="llt",.(sider.meddra_name_umls, lev="llt")]
  tmp = merge(lev.llt, sdrmdr[meddra_level=="llt",.(sider.meddra_name_umls,meddra)], by="sider.meddra_name_umls") 
  tmp = unique(merge(tmp, mdr.wide[primary_soc_fg=="Y",.(meddra=llt_code,pt_code)], by="meddra")[,.(sider.meddra_name_umls,lev="pt", meddra=pt_code)])
  gt1 = as.data.table(table(tmp$sider.meddra_name_umls))[N>1]$V1 #umls mapped to >1 pt; just pick one (we want 1:1 mapping)
  map = tmp[!sider.meddra_name_umls %in% gt1]
  for(i in gt1){
      mylist = tmp[sider.meddra_name_umls == i]$meddra
      add = fixmupt( i, mylist, mdr.wide)[,.(sider.meddra_name_umls=inputcode,lev,meddra)]
      map = rbind(map, add)
  }
  found = c(map$sider.meddra_name_umls)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for terms that map to a pt code, use that code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp.pt = lev[grep("pt", levels)]$sider.meddra_name_umls
  lev.pt = lev[sider.meddra_name_umls %in% tmp.pt,.(sider.meddra_name_umls, lev="pt")]
  tmp = merge(lev.pt, sdrmdr[meddra_level=="pt",.(sider.meddra_name_umls,meddra)], by="sider.meddra_name_umls") 
  gt1 = as.data.table(table(tmp$sider.meddra_name_umls))[N>1]$V1 #umls mapped to >1 pt; just pick one (we want 1:1 mapping)
  map = rbind(map, tmp[!sider.meddra_name_umls %in% gt1])
  for(i in gt1){
      mylist = tmp[sider.meddra_name_umls == i]$meddra
      add = fixmupt( i, mylist, mdr.wide)[,.(sider.meddra_name_umls=inputcode,lev,meddra)]
      map = rbind(map, add)
  }
  found = c(map$sider.meddra_name_umls)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for terms that map to a hlt code, use that code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  lev = lev[!sider.meddra_name_umls %in% found]
  tmp.hlt = lev[grep("hlt", levels)]$sider.meddra_name_umls
  lev.hlt = lev[sider.meddra_name_umls %in% tmp.hlt,.(sider.meddra_name_umls, lev="hlt")]
  tmp = merge(lev.hlt, sdrmdr[meddra_level=="hlt",.(sider.meddra_name_umls,meddra)], by="sider.meddra_name_umls")
  gt1 = as.data.table(table(tmp$sider.meddra_name_umls))[N>1]$V1
  map = rbind(map, tmp[!sider.meddra_name_umls %in% gt1])
  for(i in gt1){
      mylist = tmp[sider.meddra_name_umls ==i]$meddra
      add = fixmupt( i, mylist, mdr.wide, start="hlt")[,.(sider.meddra_name_umls=inputcode,lev,meddra)]
      map = rbind(map, add)
  }
  found = c(map$sider.meddra_name_umls)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for terms that map to a hlgt code, use that code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  lev = lev[!sider.meddra_name_umls %in% found]
  tmp.hlgt = lev[grep("hlgt", levels)]$sider.meddra_name_umls
  lev.hlgt = lev[sider.meddra_name_umls %in% tmp.hlgt,.(sider.meddra_name_umls, lev="hlgt")]
  tmp = merge(lev.hlgt, sdrmdr[meddra_level=="hlgt",.(sider.meddra_name_umls,meddra)], by="sider.meddra_name_umls")
  gt1 = as.data.table(table(tmp$sider.meddra_name_umls))[N>1]$V1
  map = rbind(map, tmp[!sider.meddra_name_umls %in% gt1])
  for(i in gt1){
      mylist = tmp[sider.meddra_name_umls ==i]$meddra
      add = fixmupt( i, mylist, mdr.wide, start="hlgt")[,.(sider.meddra_name_umls=inputcode,lev,meddra)]
      map = rbind(map, add)
  }
  found = c(map$sider.meddra_name_umls)
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for terms that map to a soc code, use that code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  lev = lev[!sider.meddra_name_umls %in% found]
  tmp.soc = lev[grep("soc", levels)]$sider.meddra_name_umls
  lev.soc = lev[sider.meddra_name_umls %in% tmp.soc,.(sider.meddra_name_umls, lev="soc")]
  tmp = merge(lev.soc, sdrmdr[meddra_level=="soc",.(sider.meddra_name_umls,meddra)], by="sider.meddra_name_umls")
  map = rbind(map, tmp)
  found = c(map$sider.meddra_name_umls)

  length(found) == uniqueN(sdrmdr$sider.meddra_name_umls) # [1] TRUE
  length(found) # [1] 49334
  uniqueN(sdrmdr$sider.meddra_name_umls) # [1] 49334

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add meddra codes to sider
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdr = merge(sdr, map, by="sider.meddra_name_umls")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # fwrite(sdr, row=F, quo=F, sep='\t', file="/sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/sider_awcQC.txt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUM STAT - frac drug/ae pairs with freq info
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  length(which(!is.na(sdr$sider.freqlo)))/length(sdr$sider.freqlo) # [1] 0.3925423

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUM STAT - n stitch ids
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uniqueN( gsub("CID1", "", sidA$stitchF) ) # [1] 1430 ... flat ids
  uniqueN( gsub("CID0", "", sidA$stitchS) ) # [1] 1556 ... stereo ids
  length( which(unique(gsub("CID1", "", sidA$stitchF)) %in% unique(gsub("CID0", "", sidA$stitchS)) )) # [1] 966 ... flat and stereo ids sharing pubchem id
  length( which(unique(gsub("CID0", "", sidA$stitchS)) %in% unique(gsub("CID1", "", sidA$stitchF)) )) # [1] 966 ... stereo ids with a pubchem id not shared by any flat id
  uniqueN(c( gsub("CID0", "", sidA$stitchS),gsub("CID1", "", sidA$stitchF))) # [1] 2020 ... ie 966+590+464 ... total unique pubchem ids



