# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unzip raw data
# 
# NOTE: Write in the path to the unzipped GitHub directory:
# "genetic_and_pharmacologic_overlap_in_schizophrenia-main" 
# that you downloaded from our GitHub page in line 9.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read -p "Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:"  user_home_directory
cd $user_home_directory
cd ./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/
unzip raw_data.zip

