# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/scratch/tmp
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  makemod = function (x) {trimws(gsub(".", " ", gsub("([.])\\1+", "\\1", tolower(make.names(x))), fixed=TRUE), "both")}
  library(data.table)   

  # Store user home directory and navigate to it
  user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
  setwd(user_home_directory)

  # Get the official home directory of the user
  home_dir <- getwd()

  # Navigate to the raw data zip file and unzip
  setwd(home_dir)
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load rxn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata")
  myrel = unique(rxrel[!is.na(RXCUI1) & !is.na(RXCUI2),.(RXCUI1,RXCUI2,RELA,SAB)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rxcuis with sab of "RXNORM"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  keep = unique(rxnLatest.conso[SAB=="RXNORM"]$RXCUI)
  length(keep) #[1] 205608
  keep = unique(c(keep, myrel[RXCUI2 %in% keep]$RXCUI1))
  length(keep) #[1] 205608 ... so, only RXCUI with SAB of RXNORM are connected to other RXCUI in the REL file

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# terms to keep 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  terms = rxnLatest.conso[ RXCUI %in% keep]$STR.mod

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write RXCUI to file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  write.table( data.table(c1=keep), col=F, row=F, quo=F, sep=" ", file="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/tmp")
