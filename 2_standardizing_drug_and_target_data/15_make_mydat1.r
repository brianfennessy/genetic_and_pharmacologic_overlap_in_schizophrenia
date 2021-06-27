# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/ind2mdr_clean.Rd
# 2. /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/adr_clean.tsv
# 3. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_rxn2rxn.Rd
# 4. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd
# 5. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/mysdr_final_map.tsv
# 6. /sc/arion/projects/psychgen/methods/rx/data/medi/mymedi_final_map.tsv
# 7. /sc/arion/projects/psychgen/methods/rx/data/meddra/v20.0/parsed_wide.tsv
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# need to load vigibase files (or the processed equivalents): drug, adr, ind, link
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Store user home directory and navigate to it
  user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
  setwd(user_home_directory)

  # Get the official home directory of the user
  home_dir <- getwd()

  # Navigate to the raw data zip file and unzip
  setwd(home_dir)

  source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")
  setwd("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load vigibase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  my.vb = load_vigibase(path=getwd(), date="sep_1_2017", filelist=c("ddx","drug","ind","link","out", "srce", "demo", "followup")) 
  drug = my.vb$drug
  ind = my.vb$ind
  ddx = my.vb$ddx
  link = my.vb$link
  outcome = my.vb$out
  srce = my.vb$srce
  demo = my.vb$demo
  followup = my.vb$followup
  rm (my.vb)
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/ind2mdr_clean.Rd")
  myadr = fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/adr_clean.tsv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load drug name mapping files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_rxn2rxn.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_vb2rxn.Rd")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load sider
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdrf="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/mysdr_final_map.tsv"
  mysdr = fread(sdrf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load medi - Update 6/26/21 - Try to get rid of anything MEDI related
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdif="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/medi/mymedi_final_map.tsv"
  mymdi = fread(mdif)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load meddra
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdrf="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/meddra/v20.0/parsed_wide.tsv"
  mymdr = fread(mdrf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# clean up/reformat some of the vigibase files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ##
  ## link
  ##
  link[ Dechallenge1 %in% c("Unknown", "Not_applicable"), Dechallenge1 := NA ]
  link[ Dechallenge2 %in% c("Effect_unknown", "Not_applicable"), Dechallenge2 := NA ]
  link[ Rechallenge1 %in% c("Unknown", "Not_applicable"), Rechallenge1 := NA ]
  link[ Rechallenge2 %in% c("Effect_unknown", "Not_applicable"), Rechallenge2 := NA ]
  link[ , TimeToOnsetMin := as.numeric(TimeToOnsetMin) ]
  link[ , TimeToOnsetMax := as.numeric(TimeToOnsetMax) ]
  
  ##
  ## outcome
  ##
  outcome[!Serious %in% c("Y","N"), Serious := NA ]
  outcome[, Seriousness := NULL ]
  colnames(outcome) = c("UMCReportId","adr.serious")
  
  ##
  ## srce
  ##
  srce[ Type == "Other", Type := NA ]
  
  ##
  ## demo
  ##
  demo[ AgeGroup=="Unknown", AgeGroup := NA ]
  demo[ Gender %in% c("Not_known","Not_applicable"), Gender := NA ]
  demo[ , DateDatabase := as.Date(DateDatabase) ]
  demo[ , FirstDateDatabase := as.Date(FirstDateDatabase) ]
  demo[ Type %in% c("Not_available_to_sender_(unknown)", "Other"), Type := NA ]
  colnames(demo) = c("UMCReportId", "demo.age", "demo.sex", "demo.lastdate", "demo.reporttype", "demo.region", "demo.firstdate")
  
  ##
  ## followup
  ##
  followup[, nFollowUpReports := 0 ]
  followup[UMCReportId != ReplacedUMCReportId, nFollowUpReports := 1 ]
  followup = followup[, list( nFollowUpReports = sum(nFollowUpReports) ), by = list(UMCReportId) ]
  
  ##
  ## adr
  ##
  myadr[ Outcome == "Unknown", Outcome:=NA ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make vigibase table with report|drug|ind|ae_role|ae in final rxnorm and meddra codes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Add WHO drug information to Vigibase -----------------
  ddx$mystuff$idmap[,Medicinalprod_Id:=as.integer(Medicinalprod_Id)]
  mydat = merge(ddx$mystuff$idmap, drug, by.x="Medicinalprod_Id", by.y="MedicinalProd_Id", allow.cartesian=TRUE)
  mydat = unique(mydat[,.(UMCReportId, Drug_Id, Basis, Substance_text.mod, Substance_Id,nSubstance=Number_of_ingredients)]) 

  ## Add meddra indication codes to vigibase ----------------
  ind[, Indication.mod := makemod(Indication) ]
  myind = unique(merge(ind, ind2mdr.clean, by.x="Indication.mod", by.y="whoind.mod", all.x=TRUE)[,.(Drug_Id,ind.meddra=meddra, ind.meddra.level=meddra_level)])

  ## Some sanity checking ----------------
  uniqueN(drug$UMCReportId) #[1] 15216085 unique reports in vigibase
  uniqueN(mydat$UMCReportId) #[1] 15216085 ... good
  uniqueN(drug$Drug_Id) #[1] 40091945 instances of a drug administered in vigibase
  uniqueN(mydat$Drug_Id) #[1] 40091945 ... good
  uniqueN(ind$Drug_Id) #[1] 17390716 instances of a drug administered in vigibase where indication is provided
  uniqueN(myind$Drug_Id) #[1] 17390716 ... good
  sum(unique(myind$Drug_Id) %in% unique(mydat$Drug_Id)) #[1] 17390716 ... sanity check, ok good all Drug_Id in ind table are in the drug table
  uniqueN(mydat$Substance_text.mod) #[1] 12552 ... number of unique who ingredients in vigibase
  uniqueN(mydat[!is.na(Basis),.(Drug_Id)]) #[1] 40047251 ... number of drug adminstrations in vigibase with adr role reported
  uniqueN(myadr$UMCReportId) #[1] 15194437 ... number of reports with adr; differs from 15216085 because earlier removed instances where meddra code for adr was 0

  ## Merge drug and ind tables ----------------
  nrow(mydat) == nrow(merge( mydat, myind, by="Drug_Id", all=TRUE)) #[1] TRUE .... ok, merging doesnt change number of rows, good
  mydat = merge( mydat, myind, by="Drug_Id", all=TRUE)

  ## Add rxcui codes ----------------
  nrow(mydat) == nrow(merge(mydat, vb2rxn$mapped.final, by=c("Substance_text.mod", "Substance_Id"), all.x=TRUE)) #[1] TRUE
  mydat = merge(mydat, vb2rxn$mapped.final, by=c("Substance_text.mod", "Substance_Id"), all.x=TRUE) #tmp_mydat.tsv

  ## Add adr --------------
  tmp0 = unique(myadr[,.(UMCReportId, adr.meddra=MedDRA_Id, adr.meddra.level=meddra_level,adr.meddra.n = nMedDRA_Id, adr.vbid=Adr_Id)])
  tmp1 = merge(tmp0, mydat, by="UMCReportId", allow.cartesian=TRUE)
  mydat = unique(tmp1)
  rm(tmp1)
  rm(tmp0)
  uniqueN(mydat$UMCReportId) #[1] 15194437 ... ok, so have all reports where we had an adr

  ## Add link ---------------
  colnames(link) = c("Drug_Id", "adr.vbid", "link.dechal1", "link.dechal2", "link.rechal1", "link.rechal2", "link.ttomin", "link.ttomax")
  mydat = merge( mydat, link, all.x=TRUE )

  ## Add outcome ------------------
  mydat = merge( mydat, outcome, by="UMCReportId", all.x=TRUE )

  ## Add followup ------------------
  mydat = merge( mydat, followup, by="UMCReportId", all.x=TRUE )

  ## Add source ---------------
  colnames(srce) = c("UMCReportId", "report.source" )
  mydat = merge( mydat, srce, all.x=TRUE )

  ## Add demo ---------------
  mydat = merge( mydat, demo, all.x=TRUE )

  ## write ----------------------
  mydatout = "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv"
  fwrite(mydat, row=F, quo=F, na="NA", sep="\t", file=mydatout)
