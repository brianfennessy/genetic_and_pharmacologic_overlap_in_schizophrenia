# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ ME
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. MUST READ:

#    + VigiBase file descriptions: https://charna02.u.hpc.mssm.edu/files/vigibase/VigiBase_Export_2017_Sep_1/Vigibase_export_Sep_1_2017-file_description.pdf
#    + Usefule caveats about the data: https://charna02.u.hpc.mssm.edu/files/vigibase/GuidelineUsingVigiBaseinStudies_draft.pdf

# 2. Clarification on use of the LINK file (email with UMC JAN2018) 

#    + AWC_to_UMC (23JAN2018): I notice that the ADRs are reported in terms of MedDRA codes, while the indications (i.e., "Indication" column in IND.txt) are the actual terms with no codes. In the supporting information about file types, it states that the terms in the IND.txt file are "can be decoded from MedDRA, ICD-8, ICD-9 or ICD-10." I have spent quite a bit of time trying to map these terms to MedDRA, but the mapping is incomplete. Using both exact and partial string mapping techniques against a variety of databases (e.g., UMLS), I am able to map about 80% of the terms to MedDRA. Has anyone been able to generate a more complete mapping of the indications to MedDRA by any chance? I am attempting to map the drug names in Vigibase to RxNorm. I was able to accomplish this using the ATC IDs, however I notice that the ATC IDs in Vigibase do not go up to level 5 (drug name) but rather stop at level 4. Have the Vigibase drug names been mapped to a standardized drug coding system at the level of the drug itself, rather than its class?

#    + UMC_to_AWC (24JAN2018): Mapping the indications to MedDRA will not work to 100% because not 100% of the reports are reported with MedDRA. For reports reported with ICD terminologies, mapping to MedDRA will most likely need manual work. I believe MedDRA did contain mappings to at least one ICD version in the past, not sure if they are still available though. The drugs are all coded with WHODrug, which is a widely used standard, and you should have the means to decode this within the package, a zip-file containing the necessary WHODrug files: who_ddx_dec_1_2017.zip.

#    + AWC_to_UMC (07FEB2018): Another question for you. The "Basis" column in the DRUG file, from what I gather, indicates the suspected role of the drug in the ADR. If my understanding of the LINK file is correct, this is where a given drug is explicitly linked to the ADR. However, there are instances where a drug has "Suspect" in the Basis column of the DRUG file, but is not linked to the ADR in the LINK file. For instance, for UMCReportId 1600 and Drug_Id 11750. Is this to be expected? Perhaps I am misunderstanding the different columns

#    + UMC_to_AWC (12FEB2018): There is only a record in LINK if there is information to be captured there, i.e. if information on Time to onset or Re-de challenge is missing, there would be no record in LINK. But your assumption is correct - only drugs which have  can be available in LINK.

#    + AWC_to_UMC (12FEB2018):  Does "Interacting" in basis column imply interaction between two drugs is likely cause of the ADR?

#    + UMC_to_AWC (12FEB2018):  Yes, "Interacting" implies that the interaction between two or more drugs caused the ADR. This is captured in a very heterogenous way, and  Interacting is only based on reported information, not on any assessment from UMC. There is most likely many more interactions captured in the data, that have not been marked that way in reporting.

# once IND are mapped to medra, should look at frequencies of disease classes in database to assess overall prevalence and if these are in accord with WHO morbidity statistics for regions where reports came from 
# where only 1 drug given, remove false positives using medi indications

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata
# 2. load_vigibase(path=getwd(), date="sep_1_2017", filelist=c("ddx"))
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
  setwd("/sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rxnorm 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##load("/sc/arion/projects/psychgen/methods/rx/data/rxnorm/01.03.2017/awc_parsed.Rdata")  
  load("/sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata")  
  rxnLatest.conso[, STR.mod:=makemod(STR)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ddx
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  my.vb = load_vigibase(path=getwd(), date="sep_1_2017", filelist=c("ddx")) 
  ddx = my.vb$ddx

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# initiate tracker (to be set when done as ddx$mystuff$sub2rxn)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mysub = unique(ddx$mystuff$idmap[!is.na(Substance_Id),.(Substance_text.mod, Substance_Id)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format rxn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  myrxn = unique(rxnLatest.conso[,.(RXCUI,SAB,CODE,TTY,STR.mod)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ddx sub -> rxn str (exact matches)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mapped = merge(mysub, myrxn, by.x="Substance_text.mod", by.y="STR.mod")
  unmapped = unique(mysub[!Substance_text.mod %in% mapped$Substance_text.mod]$Substance_text.mod)
  mapped$method = "exact"
  length(unmapped) #[1] 5478
  uniqueN(mapped$Substance_text.mod) #[1] 7570

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ddx mp -> rxn str (exact matches)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mymed = unique(ddx$mystuff$idmap[Substance_text.mod %in% unmapped & Number_of_ingredients==1,.(Medicinalprod_text.mod, Substance_text.mod, Substance_Id)])
  tmp = unique(merge(mymed, myrxn, by.x="Medicinalprod_text.mod", by.y="STR.mod")[,.(Substance_text.mod,Substance_Id,RXCUI,SAB,CODE,TTY)])
  tmp$method = "exact"
  mapped = rbind(mapped, tmp)
  unmapped = unique(mysub[!Substance_text.mod %in% mapped$Substance_text.mod]$Substance_text.mod)
  length(unmapped) #[1] 4874
  uniqueN(mapped$Substance_text.mod) #[1] 8174

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ddx sub -> rxn str (remove trailing salt names (eg, "sodium") and try exact matching again)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  unmapped.dt = data.table(Substance_text.mod=unmapped, Substance_text.mod.trim=rmsalts(unmapped))
  unmapped.dt[Substance_text.mod.trim=="",Substance_text.mod.trim:=NA]
  found = merge(merge(mysub, unmapped.dt, by="Substance_text.mod"), myrxn, by.x="Substance_text.mod.trim", by.y="STR.mod")
  found[, method:="exact.trimmed"]
  found = found[, colnames(mapped), with=F]
  mapped = rbind(mapped, found)
  unmapped = unique(mysub[!Substance_text.mod %in% mapped$Substance_text.mod]$Substance_text.mod)
  length(unmapped) #[1] 3812
  uniqueN(mapped$Substance_text.mod) #[1] 9236

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add who mp terms with 1 ingred to unmapped list for partial matching 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  unmapped.dt = data.table(Substance_text.mod=unmapped, Substance_text.mod.trim=rmsalts(unmapped))
  unmapped.dt[Substance_text.mod.trim=="",Substance_text.mod.trim:=NA]
  unmapped = unique(c(unmapped.dt$Substance_text.mod.trim, 
                      unique(ddx$mystuff$idmap[Substance_text.mod %in% unmapped & Number_of_ingredients==1]$Medicinalprod_text.mod)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run partial match through rxnorm api
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scratch_dir = "/sc/arion/projects/psychgen/methods/rx/scratch"
  apx = get_rxcui_approx(unmapped, scratch_dir) ## 13apr2021: throwing an error, appears to be due to ability to access internet in foreach

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge output of partial mapping with rxnorm so additional rxnorm info is viewable to evaluate matches
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  apx = unique(merge(apx, rxnLatest.conso, by.x="rxcui", by.y="RXCUI", allow.cartesian=TRUE)[,.(who=input, rxn=STR.mod, rxcui, score=score/100, SAB, TTY)])
  apx = apx[!is.na(who)]
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge output of partial mapping with who data so more who info is available to evaluate matches 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sbids = unique(apx$who)[unique(apx$who) %in% ddx$mystuff$idmap$Substance_text.mod.trim]
  mpids = unique(apx$who)[unique(apx$who) %in% ddx$mystuff$idmap$Medicinalprod_text.mod & !unique(apx$who) %in% sbids]
  tmp1 = merge(apx, unique(ddx$mystuff$idmap[,.(who.sub=Substance_text.mod,Substance_text.mod.trim,Substance_Id, who.mp=NA, text="sub")]), 
                by.x="who",  by.y="Substance_text.mod.trim")
  tmp2 = merge(apx, unique(ddx$mystuff$idmap[Medicinalprod_text.mod %in% mpids,.(Medicinalprod_text.mod,who.sub=Substance_text.mod,who.sub.trim=Substance_text.mod.trim,Substance_Id, text="mp")]), by.x="who",  by.y="Medicinalprod_text.mod")
  colnames(tmp1)[colnames(tmp1)=="who"] = "who.sub.trim"
  colnames(tmp2)[colnames(tmp2)=="who"] = "who.mp"
  apx = rbind(tmp1, tmp2, use.names=TRUE) 
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ASSESS partial match through rxnorm api : scores 50 or more
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  head(unique(apx[score>0.5,.(who.sub,rxn)]),50) # ... shows they are all correct
  uniqueN(apx[score>0.5]$who.sub) == uniqueN(apx[score>0.5 & SAB=="RXNORM"]$who.sub) #[1] TRUE 
  #
  # ... shows all good matches have an RXCUI already in RXNORM SAB ... just keep those
  #
  mapped = rbind(mapped, unique(apx[score>0.5 & SAB=="RXNORM",.(Substance_text.mod=who.sub, Substance_Id, RXCUI=rxcui, SAB, CODE=rxcui, TTY, method="rxnormapprox_scoregt50")]))
  unmapped = unique(mysub[!Substance_text.mod %in% mapped$Substance_text.mod]$Substance_text.mod)
  apx = apx[who.sub %in% unmapped]
  uniqueN(mapped[method=="rxnormapprox_scoregt50"]$Substance_text.mod) #[1] 1408
  uniqueN(mapped[ method %in% c("exact","exact.trimmed","rxnormapprox_scoregt50") ]$Substance_text.mod)/13048 #[1] 0.7916156

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# partial match through rxnorm api : scores 50 or less but with other parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  library(stringr)
  apx[ ,who.sub.firstword := tstrsplit(who.sub, split=" ", fixed=TRUE, keep=1L)] #extract first words of who and rxn strings
  apx[ ,who.sub.trim.firstword := tstrsplit(who.sub.trim, split=" ", fixed=TRUE, keep=1L)]
  apx[ ,who.mp.firstword := tstrsplit(who.mp, split=" ", fixed=TRUE, keep=1L)]
  apx[ ,rxn.firstword := tstrsplit(rxn, split=" ", fixed=TRUE, keep=1L)]
  apx[ ,who.sub.nword:=str_count(who.sub, " ")+1] #count nwords in who and rxn strings
  apx[ ,who.sub.trim.nword:=str_count(who.sub.trim, " ")+1]
  apx[ ,who.mp.nword:=str_count(who.mp, " ")+1]
  apx[ ,rxn.nword:=str_count(rxn, " ")+1]
  for(i in 1:nrow(apx)){ 
      cat("\r", i, "of", nrow(apx), "\t\t")  #count overlapping words between rxn and who strings
      str1s = unlist(strsplit(apx$who.sub[i], split=" "))
      str1t = unlist(strsplit(apx$who.sub.trim[i], split=" "))
      str1m = unlist(strsplit(apx$who.mp[i], split=" "))
      str2 = unlist(strsplit(apx$rxn[i], split=" "))
      apx[i, who.sub.olap.nword := length(intersect(str1s, str2))]
      apx[i, who.sub.trim.olap.nword := length(intersect(str1t, str2))]
      apx[i, who.mp.olap.nword := length(intersect(str1m, str2))]
  }
  apx[, who.sub.olap.pctwho:=who.sub.olap.nword/who.sub.nword] #percent of who and rxn words that are in overlap
  apx[, who.sub.trim.olap.pctwho:=who.sub.trim.olap.nword/who.sub.trim.nword]
  apx[, who.mp.olap.pctwho:=who.mp.olap.nword/who.mp.nword]
  apx[, who.sub.olap.pctrxn:=who.sub.olap.nword/rxn.nword]
  apx[, who.sub.trim.olap.pctrxn:=who.sub.trim.olap.nword/rxn.nword]
  apx[, who.mp.olap.pctrxn:=who.mp.olap.nword/rxn.nword]
  who.sub.good = apx[who.sub.olap.nword>0 & 
                     who.sub.firstword==rxn.firstword & 
                     who.sub.olap.pctwho>=0.5 & 
                     who.sub.olap.pctrxn>=0.5] #same first word and >50% of words shared
  who.sub.trim.good = apx[who.sub.trim.olap.nword>0 & 
                     who.sub.trim.firstword==rxn.firstword & 
                     who.sub.trim.olap.pctwho>=0.5 & 
                     who.sub.trim.olap.pctrxn>=0.5] #same first word and >50% of words shared
  test100.sub = sample(unique(who.sub.good$who.sub), 100)
  test100.sub.trim = sample(unique(who.sub.trim.good$who.sub), 100)
  test100.dt.sub = c()
  test100.dt.sub.trim = c()
  for (i in 1:100) {
      new2 = head(who.sub.trim.good[who.sub == test100.sub.trim[i],.(who.sub,rxn)],1)
      test100.dt.sub.trim = rbind(test100.dt.sub.trim,new2) 
  }
  mapped = rbind(mapped, unique(who.sub.trim.good[,.(Substance_text.mod=who.sub, Substance_Id, RXCUI=rxcui, SAB, CODE=rxcui, TTY, method="rxnormapprox_scorelt50_samefw_share50pctwords")]))
  unmapped = unique(mysub[!Substance_text.mod %in% mapped$Substance_text.mod]$Substance_text.mod)
  apx = apx[who.sub %in% unmapped]
  apx = merge(apx, data.table(who.sub=unmapped), all.y=TRUE) 
  uniqueN(mapped[method=="rxnormapprox_scorelt50_samefw_share50pctwords"]$Substance_text.mod) #[1] 560
  uniqueN(mapped[ method %in% c("exact","exact.trimmed","rxnormapprox_scoregt50","rxnormapprox_scorelt50_samefw_share50pctwords")]$Substance_text.mod)/13048 #[1] 0.834534

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vb2rxn = list( "all" = mysub, "mapped"=mapped, "unmapped"=unmapped )
  # save(vb2rxn, file="/sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd")
