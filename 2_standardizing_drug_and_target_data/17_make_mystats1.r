# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv
# 2. /sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/mysdr_final_map.tsv
# 3. /sc/arion/projects/psychgen/methods/rx/data/medi/mymedi_final_map.tsv
#  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat_report_stats.tsv
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

  # Navigate to the unzipped raw data directory per the user input
  setwd(home_dir)
  source("./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/myfunctions.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load vigibase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mydatout = "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv"
  mydat = fread(mydatout)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load sider
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdrf="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/sider4.1/mysdr_final_map.tsv"
  mysdr = fread(sdrf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load medi - Update 6/26/21 - Try to get rid of MEDI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdif="./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/medi/mymedi_final_map.tsv"
  mymdi = fread(mdif)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load meddra
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
  mymdr = copy(mdr.wide)
  mymdr[,llt_name:=NULL]
  mymdr[,llt_code:=NULL]
  mymdr = unique(mymdr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get stats from mydat
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mystats = mydat_stats(mydat, mysdr, mymdi, mymdr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  outf = "./genetic_and_pharmacologic_overlap_in_schizophrenia-main/1_inhouse_functions_and_setup/raw_data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat_report_stats.tsv"  
  fwrite(mystats, na="NA", outf)

