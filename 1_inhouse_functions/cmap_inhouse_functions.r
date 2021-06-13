# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. N/A
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. N/A
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# functions for analyzing cmap data released with the 2017 Cell paper on l1000 data (Subramanian et al, PMID: 29195078)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rhdf5)
library(methods)
library(data.table)
library(cmapR)
library(foreach)
library(parallel)
library(doMC)
library(corrplot)
library(Hmisc)
library(variancePartition)
options(cores = parallel:::detectCores())
registerDoMC(15)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parse cmap metadata files - input is the directory you want the parsed Rdata file to go
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_parsemeta <- function(datadir){

    filename <- paste(datadir, "cmap_parsed_metadata.Rdata", sep="/")

    if (!file.exists(filename)){

        cat("No parsed metadata file exists in specified directory\n")
        mymet <- list()

        cat("...parsing gene metadata\n")
          mymet[["rmeta"]] <- read.delim(paste(datadir,"GSE70138_Broad_LINCS_gene_info_2017-03-06.txt", sep="/"), sep="\t", row.names=1, stringsAsFactors=F)
          mymet[["rmeta"]]$gene_id <- rownames(mymet[["rmeta"]])

        cat("...parsing drug metadata\n")
          pert <- read.delim(paste(datadir,"GSE70138_Broad_LINCS_pert_info_2017-03-06.txt", sep="/"), sep="\t", row.names=1, stringsAsFactors=F)
          pert$pert_id <- rownames(pert)

        cat("...parsing cell metadata\n")
          cell <- read.delim(paste(datadir,"GSE70138_Broad_LINCS_cell_info_2017-04-28.txt", sep="/"), sep="\t", row.names=1, stringsAsFactors=F)
          cell["MNEU.E","primary_site"] <- "in vitro nervous system"
          cell["NEU.KCL","primary_site"] <- "in vitro nervous system"
          cell["NEU","primary_site"] <- "in vitro nervous system"
          cell["HUES3", "primary_site"] <- "embryo"         
          cell$cell_id <- rownames(cell)

        cat("...parsing individual experiment metadata (levels 1-4)\n")
          cmet <- read.delim(paste(datadir,"GSE70138_Broad_LINCS_inst_info_2017-03-06.txt", sep="/"), sep="\t", row.names=1, stringsAsFactors=F)
          cmet$experiment_id <- rownames(cmet)
          cmet <- merge(cmet, cell, by="cell_id")
          cmet <- merge(cmet,pert,by=c("pert_id","pert_type"))
          mymet[["cmeta.l234"]] <- cmet

        cat("...parsing signature metadata (level 5)\n")        
          sigi <- read.delim(paste(datadir,"GSE70138_Broad_LINCS_sig_info_2017-03-06.txt", sep="/"), sep="\t", stringsAsFactors=F)
          sigm <- read.delim(paste(datadir,"GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt", sep="/"), sep="\t", row.names=1, stringsAsFactors=F)
          cmet <- merge(sigi, sigm, by=c("sig_id","pert_type","pert_id"), suffixes=c(".phase2",".phase1"))
          cmet <- merge(cmet, cell, by="cell_id")
          cmet <- merge(cmet,pert,by=c("pert_id","pert_type"))
          rownames(cmet) <- cmet$sig_id
          mymet[["cmeta.lev5"]] <- cmet

        cat("...saving parsed metadata to", filename, "\n")
          save(mymet, file = filename)

    } else {

        cat("Parsed metadata file already exists in specified directory, loading...\n")
          load(filename)

    }

    cat("DONE\n")
    return(mymet)

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load cmap expression data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_loadexp <- function(datadir, level) {
    filelist <- list.files(datadir, pattern="GSE70138_Broad_LINCS_Level", full.names=T)
    names(filelist) <- paste("Level", tstrsplit(tstrsplit(filelist, split="Level", keep=2L)[[1]], split="_", keep=1L)[[1]], sep="")
    if (level==2) mydat <- parse.gctx(filelist["Level2"])
    if (level==3) mydat <- parse.gctx(filelist["Level3"])
    if (level==4) mydat <- parse.gctx(filelist["Level4"])
    if (level==5) mydat <- parse.gctx(filelist["Level5"])
    return(mydat)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# annotate expression data using metadata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_annotate <- function(mymet, mydat, level) {
    cat("Annotating rows (aka genes)...\n")
    mydat <- annotate.gct(mydat, mymet$rmeta, dim="row", keyfield="gene_id")
    cat("Annotating columns (aka data)...\n")
    if(level<5){
        mydat <- annotate.gct(mydat, mymet$cmeta.l234, dim="col", keyfield="experiment_id")
    } else if (level==5) {
        mydat <- annotate.gct(mydat, mymet$cmeta.lev5, dim="col", keyfield="sig_id")
    } else {
        mydat <- annotate.gct(mydat, mymet$cmeta.levX, dim="col", keyfield="aggregate_id")
    }
    return(mydat)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subset annotated expression data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_subset <- function(mydat, cvar=NULL, rvar=NULL){ 
    c.idx <- r.idx <- c()
    if(!is.null(cvar)) for (i in 1:length(cvar)) c.idx <- unique(c(c.idx,which(mydat@cdesc[, names(cvar)[i]] %in% cvar[[i]])))
    if(!is.null(rvar)) for (i in 1:length(rvar)) r.idx <- unique(c(r.idx,which(mydat@rdesc[, names(rvar)[i]] %in% rvar[[i]])))    
    if (!is.null(c.idx) & !is.null(r.idx)){
        mydat_sub <- subset.gct(mydat, cid=c.idx, rid=r.idx)
    } else if (!is.null(c.idx) & is.null(r.idx)) {
        mydat_sub <- subset.gct(mydat, cid=c.idx)    
    } else if (is.null(c.idx) & !is.null(r.idx)) {
        mydat_sub <- subset.gct(mydat, rid=r.idx)    
    }
    return(mydat_sub)
}
    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# map level234 identifiers to level5
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_idmap <- function(mymet){
    dt <- as.data.table(mymet$cmeta.lev5)[,.(primary_site,sig_id,distil_id)]
    dt <- dt[, list( distil_id = unlist(strsplit(distil_id, "|", fixed=T))), by=list(sig_id)]
    return(dt)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate consensus expression values by combining across experiments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_consensus <- function(mydat, ids, level) {

    ## - for general explanation of this procedure, see "Level 5 - Replicate-consensus signatures (MODZ)" section in methods of Subramanian et al Cell paper (2017)
    ## - for full details of procedure, though, it was necessary to go through the cmap Matlab source code at:
    ##    https://github.com/cmap/cmapM/blob/master/%2Bcmapm/%40Pipeline/level4_to_level5.m
    ##    https://github.com/cmap/cmapM/blob/master/%2Bcmapm/%40Pipeline/private/clip.m
    ##    https://github.com/cmap/cmapM/blob/master/%2Bcmapm/%40Pipeline/private/modzs.m
    
    ## subset expression matrix for experiment ids (keeps all genes)
    if(level==4) mydat_sub <- cmap_subset(mydat, cvar=list("experiment_id"=ids))
    if(level==5) mydat_sub <- cmap_subset(mydat, cvar=list("sig_id"=ids))
    
    ## set large outlier values to hi/lo thresholds
    mydat_sub@mat[mydat_sub@mat<=-10] <- -10
    mydat_sub@mat[mydat_sub@mat>=10] <- 10
    
    ## subset expression matrix for landmark genes and experiment ids
    if(level==4) ex <- cmap_subset(mydat, rvar=list("pr_is_lm"=1), cvar=list("experiment_id"=ids))
    if(level==5) ex <- cmap_subset(mydat, rvar=list("pr_is_lm"=1), cvar=list("sig_id"=ids))
    
    ## get correlations between experiment replicates (landmark genes only)
    ex.cor <- cor(ex@mat, method="spearman")
    
    ## set self-correlations to zero 
    diag(ex.cor) <- 0
    
    ## set negative correlations to zero 
    ex.cor[ex.cor<=0] <- 0
    
    ## weight replicates by summing its correlations with other replicates (divide by 2 since matrix is symmetric)
    wt <- colSums(ex.cor)*0.5
    
    ## set any weights less than 0.01 to 0.01
    wt[wt<=0.01] <- 0.01
    
    ## normalize weights by dividing by the sum of the weights
    sum_wt <- sum(abs(wt))
    weights <- wt/sum_wt 
    
    ## transform expression values for each replicate by multiplying them by the weight for this replicate
    for(i in names(weights)){       
        mydat_sub@mat[,i] <- mydat_sub@mat[,i]*weights[i]
    }
    
    ## generate consensus score for the experiment by summing the transformed values for its replicates
    cons <- rowSums(mydat_sub@mat)
    
    return(cons)
    
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# aggregate exp matrix by annotations (ie, tissue, dose, etc)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_aggregate <- function(mydat, mymet, level, vars=NULL, test=FALSE){

    ## check that either level 4 or 5 is given as argument since this doesnt make sense for other levels
    if (!(level %in% c(4,5))){
        stop("ERROR: level argument must be set to either 4 or 5\n")
    }

    ## aggregating level 4 -> level 5 (sanity check that we can reproduce level5 data supplied in data release)
    if (level==4){              
        idmap <- cmap_idmap(mymet) #maps level4 "experiment_id" with level5 "sig_id"
        it <- unique(idmap$sig_id)
        if(test) it <- it[1:100]
        res_final <- foreach (i = 1:length(it), .combine=cbind)%dopar%{
            if (i %% 10 == 0 ) cat("\r",i," of ", length(it),"\t\t")
            ids <- idmap[sig_id==it[i]]$distil_id
            ret <- cmap_consensus(mydat, ids, level) #this is the meat of the function, see comments above
            ret
        }
        colnames(res_final) <- it
    }

    ## aggregating level 5 -> broader groupings (ie, combining all experiments by tissue/drug)
    if (level==5){

        ## define variables that we are allowed to aggregate on 
        allowed_vars <- c("primary_site", "cell_id","sample_type", "pert_iname", "pert_idose", "pert_itime")

        ## checks
        if (sum(vars %in% allowed_vars)==0 | sum(!(vars %in% allowed_vars))>0) stop("vars supplied are not in [ ", capture.output(cat(allowed_vars))," ]\n")        
        if (is.null(vars)) stop("ERROR: level argument set to 5 but no vars argument\n")
        
        idmap <- as.data.table(mymet$cmeta.lev5) #split metadata according to the variables we are aggregating on 
        master <- split(idmap, by=vars)
        if (test) master <- master[1:100]
        master.names <- unlist(lapply(master,function(x) paste(unlist(unique(x[,vars,with=F])), collapse="|")))
        mymat <- foreach (i = 1:length(master), .combine=cbind)%dopar%{
            if (i %% 10 == 0 ) cat("\r",i," of ", length(master),"\t\t")
            ids <- master[[i]]$sig_id
            if (length(ids)==1) {
                ret <- mydat@mat[,ids]
            } else {
                ret <- cmap_consensus(mydat, ids, level)                
            }
            ret
        }
        cat("\n")
        master.names.index <- as.integer(sub("result.", "", colnames(mymat)))
        master.names <- master.names[master.names.index]
        colnames(mymat) <- master.names
        out <- cmap_makegct(mydat, mymet, mymat, vars, idmap)
    }
    return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# aggregate level 5
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_aggregate_l5_drugset <- function(mydat, mymet, drug.start.index, drug.stop.index, vars=NULL){

    ## aggregating level 5 -> broader groupings (ie, combining all experiments by tissue/drug)
    level <- 5

    ## define variables that we are allowed to aggregate on
    allowed_vars <- c("primary_site", "cell_id","sample_type", "pert_iname", "pert_idose", "pert_itime")

    ## checks
    if (sum(vars %in% allowed_vars)==0 | sum(!(vars %in% allowed_vars))>0) {
        stop("vars supplied are not in [ ", capture.output(cat(allowed_vars))," ]\n")
    }
    if (is.null(vars)) {
        stop("ERROR: level argument set to 5 but no vars argument\n")
    }
    idmap <- as.data.table(mymet$cmeta.lev5) #split metadata according to the variables we are aggregating on
    master <- split(idmap, by=vars)
    master <- master[drug.start.index:drug.stop.index]
    master.names <- unlist(lapply(master,function(x) paste(unlist(unique(x[,vars,with=F])), collapse="|")))
    mymat <- foreach (i = 1:length(master), .combine=cbind)%dopar%{
        if (i %% 10 == 0 ) cat("\r",i," of ", length(master),"\t\t")
        ids <- master[[i]]$sig_id
        if (length(ids)==1) {
            ret <- mydat@mat[,ids]
        } else {
            ret <- cmap_consensus(mydat, ids, level)
        }
        ret
    }
    cat("\n")
    master.names.index <- as.integer(sub("result.", "", colnames(mymat)))
    master.names <- master.names[master.names.index]
    colnames(mymat) <- master.names
    out <- cmap_makegct(mydat, mymet, mymat, vars, idmap)
    return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make gct object from aggregated expression matrix
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_makegct <- function(mydat, mymet, mymat, vars, idmap){
    annovars <- c("primary_site", "cell_type", "cell_id",  "sample_type", "pert_iname", "pert_id", "pert_itime", "pert_idose", "distil_id", "sig_id", "n_sig_id")
    mydat@mat <- mymat
    rows <- rownames(mymat)
    cols <- colnames(mymat)
    mydat@rid <- rows
    mydat@cid <- cols
    mydat@rdesc <- data.frame("id"=rows)
    mydat@cdesc <- data.frame("id"=cols)
    myanno <- copy(idmap)
    myanno[, aggregate_id := Reduce(function(...) paste(..., sep = "|"), .SD), .SDcols = vars]
    myanno[, n_sig_id := uniqueN(sig_id), by=aggregate_id]
    myanno <- myanno[, lapply(.SD, function(x) paste(sort(unique(x)),collapse = "|")), by=aggregate_id, .SDcols=annovars]    
    mymet$cmeta.levX <- myanno
    mydat <- cmap_annotate(mymet, mydat, level=10)
    return(mydat)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate enrichment score
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmap_es <- function(mydat, p=1, gset, direction){
    count <- 0
    mdat_rank <- rank.gct(mydat, dim="column") # rank genes by diff expr in experiment
    keepers <- c()
    for (col in colnames(mdat_rank@mat)) {
        deltas <- c()
        dt <- data.table(genes=rownames(mdat_rank@mat), rank=mdat_rank@mat[,col])[order(rank)]
        L <- merge(dt, as.data.table(mydat@rdesc), by.x="genes", by.y="id")[order(rank)]$pr_gene_symbol #rank-ordered list of genes
        x <- which(L %in% gset$gene)
        last <- x[length(x)] - 1
        p <- 1
        for (i in 1:last){
            count <- count + 1
            cat("\r", count,"of",ncol(mdat_rank@mat)*last,"\t\t")
            overlap <- l[l %in% gset$gene]
            n_r <- sum(abs(gset$z^p))
            p_hit  <- sum(abs(gset[gene %in% overlap]$z)^p/n_r)
            p_mis <- (1 / (length(l) - (uniqueN(gset$gene)-length(overlap)))) * (uniqueN(gset$gene)-length(overlap))
            deltas <- c(deltas, p_hit - p_mis)
        }
        keepers <- c(keepers,paste("i=",which(deltas==max(deltas)),";es=",max(deltas),sep=""))
        names(keepers)[length(keepers)] <- col
    }
    return(keepers)
}

