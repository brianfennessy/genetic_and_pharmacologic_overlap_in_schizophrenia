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
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/sider_awcQC.txt
# 2. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/meddra_all_indications.tsv
# 3. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/meddra.tsv
# 4. mdr = parsemdr(version="v20.0")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/sider_ind_awcQC.txt
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SIDER v4.1 uses meddra 16.1 
### (per http://sideeffects.embl.de/ - "The current version (SIDER 4.1) has been released on October 21, 2015. This release uses the MedDRA dictionary (version 16.1)")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  setwd("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects")
  source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider ae
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdr.ae = fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/sider_awcQC.txt")
  sdr.ae = unique(sdr.ae[,.(sider.meddra_name_umls,meddra,lev)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# indications
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdr.ind = unique(fread("sider4.1/meddra_all_indications.tsv", na=c("","NA"), header=F)[!is.na(V6),.(sider.stitchF=V1, sider.meddra_name_umls=V6)])
  sdr.ind.ok = sdr.ae[sider.meddra_name_umls %in% sdr.ind$sider.meddra_name_umls]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# map indications not in ae to final meddra
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdrmdr = unique(fread('sider4.1/meddra.tsv', header=F, sep='\t')[,.(sider.meddra_name_umls=V1,meddra=V3)])
  sdr.ind = merge(sdr.ind, sdrmdr)
  sdr.ind.notok = unique(sdr.ind[!sider.meddra_name_umls %in% sdr.ind.ok$sider.meddra_name_umls,.(sider.meddra_name_umls, meddra)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# meddra
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
# map meddra terms supplied by sider (meddra v16) to meddra v20 pt level (and above)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdrmdr = unique(merge(mdr, sdr.ind.notok, by="meddra")[,.(sider.meddra_name_umls,meddra,meddra_level)])
  lev = merge(mdr, sdr.ind.notok, by="meddra")[, list( levels = paste(unique(meddra_level), collapse="|")), by=list(sider.meddra_name_umls) ]
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for terms that are supplied as llt level only, map to pt
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  lev.llt = lev[levels=="llt",.(sider.meddra_name_umls, lev="llt")]
  tmp = merge(lev.llt, sdrmdr[meddra_level=="llt",.(sider.meddra_name_umls,meddra)], by="sider.meddra_name_umls") 
  tmp = unique(merge(tmp, mdr.wide[primary_soc_fg=="Y",.(meddra=llt_code,pt_code)], by="meddra")[,.(sider.meddra_name_umls,lev="pt", meddra=pt_code)])
  gt1 = as.data.table(table(tmp$sider.meddra_name_umls))[N>1]$V1 #umls mapped to >1 pt; just pick one (we want 1:1 mapping)
  map = tmp[!sider.meddra_name_umls %in% gt1]
  add = c()
  for(i in gt1){
      mylist = tmp[sider.meddra_name_umls == i]$meddra
      add = rbind(add, fixmupt( i, mylist, mdr.wide) )
  }
  map = rbind(map, add[,.(sider.meddra_name_umls=inputcode, lev, meddra)])
  found = c(map$sider.meddra_name_umls)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for terms that map to a pt code, use that code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp.pt = lev[grep("pt", levels)]$sider.meddra_name_umls
  lev.pt = lev[sider.meddra_name_umls %in% tmp.pt,.(sider.meddra_name_umls, lev="pt")]
  tmp = merge(lev.pt, sdrmdr[meddra_level=="pt",.(sider.meddra_name_umls,meddra)], by="sider.meddra_name_umls") 
  gt1 = as.data.table(table(tmp$sider.meddra_name_umls))[N>1]$V1 #umls mapped to >1 pt; just pick one (we want 1:1 mapping)
  map = rbind(map, tmp[!sider.meddra_name_umls %in% gt1])
  add = c()
  for(i in gt1){
      mylist = tmp[sider.meddra_name_umls == i]$meddra
      add = rbind(add, fixmupt( i, mylist, mdr.wide) )
  }
  map = rbind(map, add[,.(sider.meddra_name_umls=inputcode, lev, meddra)])
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
  found = c(map$sider.meddra_name_umls)
  length(found) == uniqueN(sdrmdr$sider.meddra_name_umls) #[1] TRUE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add meddra codes to sider
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdr.ind.ok = rbind(sdr.ind.ok, map)
  sdr.ind = unique(sdr.ind[,.(sider.meddra_name_umls, sider.stitchF)])
  sdr.ind = merge(sdr.ind, sdr.ind.ok, by="sider.meddra_name_umls")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fwrite(sdr.ind, row=F, quo=F, sep='\t', file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/sider_ind_awcQC.txt")


