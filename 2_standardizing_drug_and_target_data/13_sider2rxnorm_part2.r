# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/sider_awcQC.txt
# 2. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/sider_ind_awcQC.txt
# 3. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_sdr2rxn.Rd
# 4. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_rxn2rxn.Rd
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/mysdr.Rdata
# 2. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/mysdr_final_map.tsv
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
  setwd("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider -- need drug | phe_type (indication or ae) | phe 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdr = fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/sider_awcQC.txt",  header=T, sep='\t')
  sdr.ind = fread("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/sider_ind_awcQC.txt", header=T, sep='\t')
  sdr[, sider.meddra_name.mod:=makemod(sider.meddra_name)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sider to rxn maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_sdr2rxn.Rd")
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/drug_id_map_rxn2rxn.Rd")
  ### ... note, because used get_rxcui_in during sider mapping we lost the multingredients when we sent those back as input when making rxn2rxn

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merging sider with drug id maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp1 = unique(sdr2rxn[,.(sider.stitchF,rxcui_input=rxcui_output)]) #stitch<>rxn
  tmp2 = rxn2rxn$ing #rxn<>rxn_ing
  tmp3 = unique(merge(tmp1, tmp2)[,.(sider.stitchF=sider.stitchF,rxcui=rxcui_output)]) #stitch<>rxn_ing

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fix entries with multiple IBD matches (these arent multi-ingredients; spot check to convince yourself)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gt1 = as.data.table(table(tmp3$sider.stitchF))[N>1]$V1
  tmp4 = tmp3[!sider.stitchF %in% gt1] 
  for (i in gt1){tmp4 = rbind(tmp4, data.table(sider.stitchF=i, rxcui = sample(tmp3[sider.stitchF==i]$rxcui,1)))}
  sdr = merge(sdr, tmp4, by="sider.stitchF") 
  sdr.cln = unique(sdr[,.(rxcui, meddra, meddra_level=lev, sider.freqlo, sider.freqhi)])
  sdr.cln = sdr.cln[,list( sider.freqlo = min(sider.freqlo,na.rm=TRUE), sider.freqhi=max(sider.freqhi,na.rm=TRUE) ), by=list(rxcui,meddra,meddra_level)]
  sdr.cln[is.infinite(sider.freqlo), sider.freqlo:=NA]
  sdr.cln[is.infinite(sider.freqhi), sider.freqhi:=NA]
  sdr.cln.ln = unique(sdr.cln[,.(rxcui,meddra,meddra_level)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merging sider with drug id maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp1 = unique(sdr[,.(sider.stitchF,sider.stitchS,rxcui_input=rxcui)]) #stitch<>rxn
  tmp2 = rxn2rxn[[1]] #rxn<>rxn_ing
  tmp3 = unique(merge(tmp1, tmp2, by="rxcui_input")[,.(sider.stitchF=sider.stitchF,sider.stitchS=sider.stitchS,rxcui=rxcui_output)]) #stitch<>rxn_ing
  sdr.ind = merge(sdr.ind, tmp4, by="sider.stitchF")
  sdr.ind.cln.ln = unique(sdr.ind[,.(rxcui,meddra,meddra_level=lev)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mysdr = list("ae.clean" = sdr.cln.ln, "ae.clean_freq"=sdr.cln, "ae.full" = sdr, "ind.full"=sdr.ind, "ind.clean"=sdr.ind.cln.ln)
  save(mysdr, file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/mysdr.Rdata")

  mytab = rbind(mysdr$ae.clean[,.(rxcui,meddra,meddra_level,meddra_type="ae")], mysdr$ind.clean[,.(rxcui,meddra,meddra_level,meddra_type="ind")])
  write.table(mytab, row=F, quo=F, sep="\t", file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/mysdr_final_map.tsv")
