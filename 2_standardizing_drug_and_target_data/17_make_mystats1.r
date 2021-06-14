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
  source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
  setwd("/sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load vigibase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mydatout = "/sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv"
  mydat = fread(mydatout)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load sider
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sdrf="/sc/arion/projects/psychgen/methods/rx/data/side_effects/sider4.1/mysdr_final_map.tsv"
  mysdr = fread(sdrf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load medi
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdif="/sc/arion/projects/psychgen/methods/rx/data/medi/mymedi_final_map.tsv"
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
  outf = "/sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat_report_stats.tsv"  
  # fwrite(mystats, na="NA", outf)

