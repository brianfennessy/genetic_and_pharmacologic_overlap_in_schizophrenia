# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. rxnorm data directory must contain {RXNCONSO,RXNSTY,RXNREL}.rrf - if not using 1/3/2017 rxnorm release, may have to edit the \"rxx\" variable in this script because it uses a helper file for renaming columns that was made using 1/3/2017 release
# 2. /sc/arion/projects/psychgen/methods/rx/data/rxnorm/rxnorm.colnames
# 3. path <- commandArgs(trailingOnly=TRUE)[1]
# 4. {path}/RXNREL.RRF
# 5. {path}/RXNCONSO.RRF
# 6. {path}/RXNSTY.RRF
# 7. https://rxnav.nlm.nih.gov/REST/rxcui/{all_rxcui[i]}/status
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. {path}/awc_parsed.Rdata
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#!/usr/bin/Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# usage information
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( length(commandArgs(trailingOnly=TRUE)) == 0 ) {
        cat("\n\n---------------------------------------------------------------------------------------------------\n\n")
        cat("usage: Rscript parse_rxnorm.r path_to_rxnorm_rrf_data_directory\n\n")
        cat("\tnote1: rxnorm data directory must contain {RXNCONSO,RXNSTY,RXNREL}.rrf\n")
        cat("\tnote2: if not using 1/3/2017 rxnorm release, may have to edit the \"rxx\" variable in this script because it uses a helper file for renaming columns that was made using 1/3/2017 release\n\n")    
        cat("---------------------------------------------------------------------------------------------------\n\n")
    }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                      
# setup 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    options(stringsAsFactors=F)
    library (data.table)
    library(XML)
    library(xml2)
    library(rvest)
    library(foreach)
    library(parallel)
    library(doMC)
    options(cores = parallel:::detectCores())
    registerDoMC(15)
    setwd("/sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rxnorm data files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rxx <-  fread("/sc/arion/projects/psychgen/methods/rx/data/rxnorm/rxnorm.colnames" , header=F )
    # path <- commandArgs(trailingOnly=TRUE)[1]
    path <- "/sc/arion/projects/psychgen/methods/rx/data/rxnorm/09.04.2018/rrf"
    rxrel <- fread(paste(path, "RXNREL.RRF", sep="/"), header=F, sep="|" )
    rxnLatest.conso <- fread(paste(path, "RXNCONSO.RRF", sep="/"), header=F, sep="|")
    rxnLatest.sty <- fread(paste(path, "RXNSTY.RRF", sep="/"), header=F, sep="|")
    rxnLatest.conso [,ncol(rxnLatest.conso)-1] <- NULL
    rxnLatest.sty   [,ncol(rxnLatest.sty  )-1] <- NULL
    colnames(rxnLatest.conso) <- strsplit( rxx[ V1=="RXNCONSO", ]$V2 , "|" , fixed=T)[[1] ]
    colnames(rxnLatest.sty  ) <- strsplit( rxx[ V1=="RXNSTY"  , ]$V2 , "|" , fixed=T)[[1] ]
    rxnLatest.conso[,STRlc := tolower(STR)]
    rxrel [,ncol(rxrel)-1] <- NULL
    colnames(rxrel) <- strsplit( rxx[rxx$V1=="RXNREL"]$V2, "|" , fixed=T)[[1]]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get status for all rxcui in data files for rxnorm (this info is not in data files themself)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    all_rxcui <- unique(rxnLatest.conso$RXCUI)
    count <- 1
    res <- c()
    res_final <- foreach ( i = 1:length(all_rxcui),combine=cbind)%dopar%{
        cat("\r", i, "of", length(all_rxcui), "\t\t")
        url <- paste("https://rxnav.nlm.nih.gov/REST/rxcui/",all_rxcui[i],"/status",sep="")
        destfile <- paste(path, "/xml/", all_rxcui[i],".xml", sep="")
        download.file(url, destfile = destfile, quiet=TRUE)
        x <- xmlParse(read_xml(destfile))
        file.remove(destfile)
        ns <- getNodeSet(x,'/rxnormdata', fun = xmlToList)
        c1 <- all_rxcui[i]
        c2 <- ns[[1]]$rxcuiStatus$status
        c3 <- c4 <- c5 <- c6 <- NA
        if ( length(ns[[1]]$rxcuiStatus) > 2 ) {
            c3 <- ns[[1]]$rxcuiStatus$minConceptGroup$minConcept$rxcui
            c4 <- ns[[1]]$rxcuiStatus$minConceptGroup$minConcept$name
            c5 <- ns[[1]]$rxcuiStatus$minConceptGroup$minConcept$tty
            if (c2 == "Remapped" ) {
                c6 <- ns[[1]]$rxcuiStatus$remappedDate}
        }
        m <- as.data.table(matrix(c(c1,c2,c3,c4,c5,c6),nrow=1,ncol=6))
        colnames(m) <- c("mapped_rxcui","mapped_rxcui_status","remapped_rxcui","remapped_name","remapped_tty","remapped_date")
        res <- rbind(res,m)
        count <- count+1
        res
    }
    all_status <- do.call(rbind,res_final)
    all_status[, mapped_rxcui:=as.integer(mapped_rxcui)]
    cat("Sent for status on ", as.character(uniqueN(all_rxcui)), " RXCUIs\n", sep="")
    cat("Received status on ", as.character(nrow(all_status)), " RXCUIs\n", sep="")
    actives <- all_status[!is.na(remapped_rxcui)]$remapped_rxcui

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get rels for all active rxcui
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    active_rels <- unique(rbind(rxrel[ RXCUI1 %in% actives & !is.na(RXCUI2),.(RXCUI1,RXCUI2,RELA) ],rxrel[ RXCUI2 %in% actives & !is.na(RXCUI1),.(RXCUI1,RXCUI2,RELA) ]))
    active_rels <- merge(active_rels, all_status[,.(mapped_rxcui,mapped_rxcui_status)], by.x="RXCUI1", by.y="mapped_rxcui")[,.(RXCUI1,RXCUI2,RELA,RXCUI1_STATUS=mapped_rxcui_status)]
    active_rels <- merge(active_rels, all_status[,.(mapped_rxcui,mapped_rxcui_status)], by.x="RXCUI2", by.y="mapped_rxcui")[,.(RXCUI1,RXCUI2,RELA,RXCUI1_STATUS,RXCUI2_STATUS=mapped_rxcui_status)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
# save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # save( list = c("rxx", "rxrel", "rxnLatest.conso", "rxnLatest.sty", "all_status", "actives", "active_rels"), file = paste(path,"awc_parsed.Rdata",sep="/") )



