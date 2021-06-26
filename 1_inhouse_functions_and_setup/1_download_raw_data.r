# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load necessary libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unzip raw data
# 
# NOTE: Write in the path to the unzipped GitHub directory:
# "genetic_and_pharmacologic_overlap_in_schizophrenia-main" 
# that you downloaded from our GitHub page in line 13.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("_________")
home_path <- getwd()
setwd(paste(home_path,"genetic_and_pharmacologic_overlap_in_schizophrenia-main","1_inhouse_functions_and_setup",sep="/"))
unzip(zipfile = "raw_data.zip")