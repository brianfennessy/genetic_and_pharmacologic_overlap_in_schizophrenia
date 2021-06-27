# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd
# 2. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_sdr2rxn.Rd
# 3. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_db2rxn.Rd
# 4. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_sc2rxn.Rd
# 5. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_cmap2rxn.Rd
# 6. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_medi2rxn.Rd
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_rxn2rxn.Rd
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

  source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in mappings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_vb2rxn.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_sdr2rxn.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_db2rxn.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_sc2rxn.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_cmap2rxn.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_medi2rxn.Rd")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine into table of rxids and source
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  map = unique(vb2rxn$mapped[,.(rxcui=RXCUI,database="vigibase")])
  map = rbind(map, unique(sdr2rxn[,.(rxcui=rxcui_output,database="sider")]))
  map = rbind(map, unique(db2rxn$db.mapped[,.(rxcui,database="drugbank")]))
  map = rbind(map, unique(sc2rxn[,.(rxcui,database="seachange")]))
  map = rbind(map, unique(cmap2rxn[,.(rxcui=RXCUI, database="cmap")]))
  map = rbind(map, unique(medi2rxn[,.(rxcui, database="medi")]))
  scratchdir = "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/"
  rxcui_list =  unique(map$rxcui)
  my.ing = get_rxcui_in(rxcui_list, scratchdir, use.foreach=F)
  my.sts = get_rxcui_status(rxcui_list, scratchdir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# what are the aliens in my.sts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  aliens = unique(my.sts[ mapped_rxcui_status=="Alien" ]$mapped_rxcui)
  as.data.table(table(map[rxcui %in% aliens]$database))[order(N)]
  ##         V1    N
  ##1:     medi    4
  ##3: vigibase 2226
  ##
  table(vb2rxn$mapped[ RXCUI %in% aliens ]$method)
  ##      exact exact.trimmed
  ##       3770           631
  ##
  ## ... some (~500) of these can be rescued by running through mapping pipeline after exact matching ... skipping for now

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rxn2rxn = list("ing"=my.ing, "status"=my.sts) 
  save( rxn2rxn, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_rxn2rxn.Rd")
 
