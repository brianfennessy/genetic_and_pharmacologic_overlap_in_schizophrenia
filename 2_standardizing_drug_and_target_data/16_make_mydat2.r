# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv
# 2. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd
# 3. /sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata
# 4. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_rxn2rxn.Rd
# 5. /sc/arion/projects/psychgen/methods/rx/files/fw_subset_vb_highvolume.Rd
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. /sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd
# 2. /sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# common drugs like "insulin" and "influenza vaccine" remain unmapped in vigibase and as a result lots of reports would be filtered. these should be mappable, so need to do manually
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/myfunctions.r")
  setwd("/sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1")
  scratchdir = "/sc/arion/projects/psychgen/methods/rx/scratch"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load vigibase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mydatin1 = "/sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv"
  mydat = fread(mydatin1)
  load("/sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load rxnorm 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("/sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/awc_parsed.Rdata")  
  rxnLatest.conso[, STR.mod:=makemod(STR)]
  myrxn = unique(rxnLatest.conso[SAB=="RXNORM",.(RXCUI, TTY, STR.mod)])
  load("/sc/arion/projects/psychgen/methods/rx/files/drug_id_map_rxn2rxn.Rd")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of reports that would currently be removed if only keeping reports where all drugs mapped
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uniqueN(mydat[ is.na(rxcui_output) & !is.na(Substance_text.mod) ]$UMCReportId) #[1] 1296786

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# unmapped terms and how many reports they appear in
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  unmapped = mydat[is.na(rxcui_output) & !is.na(Substance_text.mod), .N , Substance_text.mod][order(N)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mapped terms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mapped = unique(mydat[!is.na(rxcui_output), .(Substance_text.mod, rxcui_output)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write to file for fuzzywuzzy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # fwrite( unmapped[,.(Substance_text.mod)], row=F, quo=F, col=F, sep="\t", file="/sc/arion/projects/psychgen/methods/rx/scratch/unmapped_vigibase.tsv" ) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# manual review of all unmapped terms appearing in >=1000 reports was performed, read in results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("/sc/arion/projects/psychgen/methods/rx/files/fw_subset_vb_highvolume.Rd")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add manual mapping to vb2rxn object and save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vb2rxn$mapped.manually = fw_subsets

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# approved matches after manual review
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  approved = fw_subsets[approve==1,.(Substance_text.mod=vigi, STR.mod=rxn)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# map approved matches to rxcui
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  approved = merge(approved, myrxn)
  approved = merge(unique(approved[,.(Substance_text.mod,RXCUI)]), myrxn)
  tmp1 = approved[TTY=="IN",.(Substance_text.mod,RXCUI,TTY)]
  tmp2 = approved[TTY!="IN" & !Substance_text.mod%in%tmp1$Substance_text.mod,.(Substance_text.mod,RXCUI,TTY)]
  approved = rbind(tmp1, tmp2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# show the vigibase terms (n=183) are not 1:1 mapping to rxn terms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uniqueN(approved[,.(Substance_text.mod,RXCUI)]$Substance_text.mod) #[1] 183 ... number of vigibase terms with an approved match and appearing in >1000 reports
  uniqueN(approved[,.(Substance_text.mod,RXCUI)]) #[1] 215 ... number of rxn terms these 183 terms map to 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make mappings 1:1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gt1 = unique(approved[,.(Substance_text.mod,RXCUI)])[, .N, Substance_text.mod][N>1]$Substance_text.mod
  eq1.in = approved [ TTY=="IN" & !Substance_text.mod %in% gt1 ]
  eq1.notin = approved [ TTY!="IN" & !Substance_text.mod %in% gt1 ]
  eq1.notin = unique(eq1.notin[,.(Substance_text.mod,RXCUI)])
  eq1.in = rbind(eq1.in, merge( eq1.notin, get_rxcui_in(unique(eq1.notin$RXCUI), scratchdir), by.x="RXCUI", by.y="rxcui_input")[,.(Substance_text.mod, RXCUI=rxcui_output, TTY="IN")]) 
  eq1.notin = eq1.notin [ !Substance_text.mod %in% eq1.in$Substance_text.mod ]
  for (i in gt1) {

    ## the approved matches for this vigibase term
    matches = approved[Substance_text.mod==i]

    ## if only 1 of the matches is an ingredient, keep it and move on
    if ( nrow(matches[TTY=="IN"]) == 1 ) eq1.in = rbind(eq1.in, matches[TTY=="IN"])

    ## if more than 1 of matches in an ingredient but only 1 of them is in the master mappings to rxnorm, keep it and move on
    if ( nrow(matches[TTY=="IN" & RXCUI %in% mapped$rxcui_output]) == 1 ) eq1.in = rbind(eq1.in, matches[TTY=="IN" & RXCUI %in% mapped$rxcui_output])

    ## if more than 1 of matches in an ingredient and more than 1 of them is in the master mappings to rxnorm, pick 1    
    if ( nrow(matches[TTY=="IN" & RXCUI %in% mapped$rxcui_output]) > 1 ) eq1.in = rbind(eq1.in, matches[TTY=="IN" & RXCUI %in% mapped$rxcui_output][1])

    ## if none of the matches are ingredients, get ingredients and choose 1 match
    if ( nrow(matches[TTY=="IN"]) == 0 ) {
        newone = merge(matches, get_rxcui_in(unique(matches$RXCUI), scratchdir)[1][,.(RXCUI=rxcui_input, rxcui_output)])[,.(Substance_text.mod,RXCUI=rxcui_output, TTY="IN")]
        eq1.in = rbind(eq1.in, newone )
    }

    ## redefine set with greater than 1 match
    gt1 = gt1[ !gt1 %in% eq1.in$Substance_text.mod ]

  }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how many reports did we save?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uniqueN(mydat[ is.na(rxcui_output) & !is.na(Substance_text.mod) & !Substance_text.mod %in% eq1.in$Substance_text.mod ]$UMCReportId) #[1] 619072 ... 1296786 - 619072 = 677714 reports saved

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add new mapping top vb2rxn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eq1.in = unique(merge(vb2rxn$all, eq1.in)[,.(Substance_text.mod,Substance_Id,rxcui_output=RXCUI)])
  vb2rxn$mapped.final = unique(rbind( vb2rxn$mapped.final, eq1.in ))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save vb2rxn 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # save ( vb2rxn, file="/sc/arion/projects/psychgen/methods/rx/files/drug_id_map_vb2rxn.Rd" )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# update and save mydat
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tmp0 = mydat[ ! Substance_text.mod %in% eq1.in$Substance_text.mod ]
  tmp1 = mydat[ Substance_text.mod %in% eq1.in$Substance_text.mod ]
  tmp1[,rxcui_output:=NULL]
  eq1.in[,Substance_Id:=as.integer(Substance_Id)]
  tmp1 = merge(tmp1, eq1.in, by=c("Substance_text.mod","Substance_Id"))
  mydat = rbind( tmp0, tmp1 ) 
  mydatout1 = "/sc/arion/projects/psychgen/methods/rx/data/side_effects/vigibase/VigiBase_Extract_Case_Level_2017_Sep_1/mydat.tsv"
  # fwrite(mydat, row=F, quo=F, na="NA", sep="\t", file=mydatout1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sanity check we have expected number of reports with missing rxcui
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uniqueN(mydat[ is.na(rxcui_output) & !is.na(Substance_text.mod)]$UMCReportId) #[1] 619072 ... ok, it is what it is


  
