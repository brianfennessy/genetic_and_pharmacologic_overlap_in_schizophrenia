# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/ind2mdr_clean.Rd
# 2. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_rxn2rxn.Rd
# 3. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd
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
  
  # Store user home directory and navigate to it
  user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
  setwd(user_home_directory)

  # Get the official home directory of the user
  home_dir <- getwd()

  # Navigate to the raw data zip file and unzip
  setwd(home_dir)

  setwd("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1")
  source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load rxnorm and meddra mapping files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/ind2mdr_clean.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_rxn2rxn.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_vb2rxn.Rd")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# resolve vigibase drugs that mapped to >1 rxcui ING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp1 = unique(vb2rxn$mapped[,.(Substance_Id, Substance_text.mod, rxcui_input=as.character(RXCUI))])
  tmp2 = rxn2rxn[[1]]
  tmp3 = unique(merge(tmp1, tmp2)[,.(Substance_text.mod, Substance_Id, rxcui_output)])
  tmp4 = unique(merge(tmp1, tmp2)[,.(Substance_text.mod, Substance_Id, rxcui_input, rxcui_output)])
  badguys = as.data.table(table(tmp3$Substance_text.mod))[N>1]$V1 
  good = tmp3[!Substance_text.mod %in% badguys ]
  bad = tmp3[Substance_text.mod %in% badguys ]
  tmp5 = unique(vb2rxn$mapped[ Substance_text.mod %in% badguys,.(Substance_text.mod,method) ])
  pickone = unique(tmp5[, .N , by=Substance_text.mod ][N==1]$Substance_text.mod)
  for(i in pickone){
      add = bad[Substance_text.mod==i]
      add = add[ sample(1:nrow(add), 1) ]
      good = rbind(good, add)
      bad = bad[ Substance_text.mod != i ]
      badguys = badguys[ badguys != i ]
  }
  ugh = unique(tmp5[, .N , by=Substance_text.mod ][N>1]$Substance_text.mod)
  tmp6 = unique(vb2rxn$mapped[Substance_text.mod %in% ugh,.(Substance_text.mod, rxcui_input=as.character(RXCUI), rxcui_input.sab=SAB, rxcui_input.tty=TTY, method)])
  tmp6 = merge(tmp4, tmp6, by=c("Substance_text.mod", "rxcui_input"), all.y=T) 
  tmp7 = unique(tmp6[ rxcui_input == rxcui_output & method == "exact",.(Substance_text.mod,rxcui_output) ])
  add = tmp7[ Substance_text.mod %in% tmp7[, .N, by=Substance_text.mod][N==1]$Substance_text.mod ]
  badguys = badguys[ !badguys %in% add$Substance_text.mod ]
  good = rbind( good, merge(tmp3, add) ) 
  bad = bad[ ! Substance_text.mod %in% good$Substance_text.mod ]
  tmp6 = unique(vb2rxn$mapped[Substance_text.mod %in% badguys,.(Substance_text.mod, rxcui_input=as.character(RXCUI), rxcui_input.sab=SAB, rxcui_input.tty=TTY, method)])
  tmp6 = merge(tmp4, tmp6, by=c("Substance_text.mod", "rxcui_input"), all.y=T) 
  tmp7 = unique(tmp6[ rxcui_input == rxcui_output & method == "exact",.(Substance_text.mod,rxcui_output) ])
  pickone = unique(tmp7$Substance_text.mod)
  for(i in pickone){
      add = tmp7[Substance_text.mod==i]
      add = add[ sample(1:nrow(add), 1) ]
      add = merge(add, tmp3)
      good = rbind(good, add)
      bad = bad[ Substance_text.mod != i ]
      badguys = badguys[ badguys != i ]
  }
  tmp6 = unique(vb2rxn$mapped[Substance_text.mod %in% badguys,.(Substance_text.mod, rxcui_input=as.character(RXCUI), rxcui_input.sab=SAB, rxcui_input.tty=TTY, method)])
  tmp6 = merge(tmp4, tmp6, by=c("Substance_text.mod", "rxcui_input"), all.y=T) 
  add = unique(tmp6[method=="exact" & !is.na(rxcui_output),.(Substance_text.mod, Substance_Id, rxcui_output)])
  badguys = badguys[ !badguys %in% add$Substance_text.mod ]
  tmp6 = unique(vb2rxn$mapped[Substance_text.mod %in% badguys,.(Substance_text.mod, rxcui_input=as.character(RXCUI), rxcui_input.sab=SAB, rxcui_input.tty=TTY, method)])
  tmp6 = merge(tmp4, tmp6, by=c("Substance_text.mod", "rxcui_input"), all.y=T) 
  tmp6 = tmp6[!is.na(rxcui_output)]
  pickone = unique(tmp6$Substance_text.mod)
  for(i in pickone){
      add = tmp6[Substance_text.mod==i]
      add = add[ sample(1:nrow(add), 1),.(Substance_text.mod, Substance_Id, rxcui_output) ]
      good = rbind(good, add)
      bad = bad[ Substance_text.mod != i ]
      badguys = badguys[ badguys != i ]
  }
  vb2rxn$mapped.final = good
  save ( vb2rxn, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_vb2rxn.Rd" )
