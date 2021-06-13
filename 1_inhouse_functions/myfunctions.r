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

suppressMessages(library(data.table))
suppressMessages(library(readxl))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(ggrepel))
suppressMessages(library(XML) )
suppressMessages(library(xml2) )
suppressMessages(library(rvest))
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(doMC))
suppressMessages(options(cores = parallel:::detectCores()))
suppressMessages(registerDoMC(32))
suppressMessages(options(stringsAsFactors=F))
cat("\nTHE FOLLOWING LIBRARIES WERE LOADED:\ndata.table\nreadxl\nggplot2\nggthemes\nggrepel\nXML\nxml2\nrvest\nforeach\nparallel\ndoMC\n\n")
cat("\nTHE FOLLOWING SETTINGS WERE DEFINED:\noptions(cores = parallel:::detectCores())\nregisterDoMC(32)\noptions(stringsAsFactors=F)\n\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Optional: source cmap functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# source("/sc/arion/projects/psychgen/methods/rx/psychosis_paper_2021/cmap_inhouse_functions.r")

main2subsid = function(info, dt, dtname, helpfiles){
    data = info[filegroup=="main" & file==dtname]
    cat("File", dtname, "has", nrow(data), "columns to parse...\n")
    cat("Of these,", nrow(data[!is.na(link)]), "will be linked to a subsidiary file...\n")
    for (i in 1:nrow(data)){
        cat("Parsing column", i,"\n")
        if (data[i]$columntype=="Integer"){
            dt[ , (data[i]$column) := as.integer(substr(V1, start = data[i]$start, stop = data[i]$stop))]
        } else if (data[i]$columntype=="Date"){
            dt[ , (data[i]$column) := as.Date(as.character(substr(V1, start = data[i]$start, stop = data[i]$stop)),"%Y%m%d")]
        } else {
            dt[ , (data[i]$column) := as.character(substr(V1, start = data[i]$start, stop = data[i]$stop))]
        }
        if (!is.na(data[i]$link)){
            cat("\tlinking to subsidiary file\n")
            subsid = fread(helpfiles[data[i]$link], sep="\t", header=F)
            subdata = info[filegroup=="subsidiary" & file==data[i]$link]
            for (j in 1:nrow(subdata)){
                subsid[ , (subdata[j]$column) := gsub(" ", "_", as.character(substr(V1, start = subdata[j]$start, stop = subdata[j]$stop)))]
            }
            subsid[, V1 := NULL]
            dt = merge(dt, subsid, by.x=data[i]$column, by.y="Code", all.x=T)
            dt[, (data[i]$column) := Text]
            dt[, Text := NULL]
        }
    }
    dt[, V1 := NULL]          
    setcolorder(dt,data$column)
    return(dt)
}

parsemdr = function(version) {
    dir = paste("/sc/arion/projects/psychgen/methods/rx/data/meddra/",version,"/MedAscii/", sep="")
    llt = fread(paste(dir,'llt.asc',sep=""), header=F , sep="$") [,c(1,2,3), with=F]
    colnames(llt) = c("llt_code","llt_name","pt_code")
    mr  = fread (paste(dir,'mdhier.asc',sep=""), col.names=c("pt_code","hlt_code","hlgt_code","soc_code","pt_name","hlt_name","hlgt_name","soc_name","soc_abbrev","null_field","pt_soc_code","primary_soc_fg" ,"null") , sep="$") 
    master = merge ( llt, mr, by="pt_code" , allow.cartesian=TRUE)
    master = master [ , c( "llt_code","pt_code", "hlt_code", "hlgt_code", "soc_code", "llt_name", "pt_name", "hlt_name", "hlgt_name", "soc_name", "soc_abbrev", "pt_soc_code", "primary_soc_fg"), with=F]
    return(master)
}

strSort = function (x) {
    sapply(lapply(strsplit(x, split=" "), sort), paste, collapse = " ")
}

fracmap = function(map) {
    nSUB.any = uniqueN(map[!is.na(Substance_text.mod.index)]$Substance_text.mod.index)
    nSUB.map = uniqueN(map[!is.na(Substance_text.mod.index) & !is.na(rxcui)]$Substance_text.mod.index)
    nMP.any = uniqueN(map[!is.na(Medicinalprod_text.mod.index)]$Medicinalprod_text.mod.index)
    nMP.map = uniqueN(map[!is.na(Medicinalprod_text.mod.index) & !is.na(rxcui)]$Medicinalprod_text.mod.index)
    fSUB = nSUB.map/nSUB.any
    fMP = nMP.map/nMP.any
    cat("## SUB mapped =", nSUB.map, "of", nSUB.any, "=",fSUB, "\nMP mapped =", nMP.map, "of", nMP.any, "=",fMP,"\n")
    cat("## Breakdown of mapping methods:\n")
    print(table(map$method))
}

updatemap = function(map, updater, meth) {
    ret = merge( map, updater[,.(term,rxcui)], by="term", all.x=TRUE )
    ret[is.na(rxcui.x) & !is.na(rxcui.y), rxcui.x:=rxcui.y] 
    ret[!is.na(rxcui.x) & is.na(method), method:=meth] 
    ret[, rxcui:=rxcui.x] 
    ret[, rxcui.y:=NULL] 
    ret[, rxcui.x:=NULL] 
    return(unique(ret))
}

makemod = function (x) {
    trimws(gsub(".", " ", gsub("([.])\\1+", "\\1", tolower(make.names(x))), fixed=TRUE), "both")
}

get_rxcui_approx = function(myrxcui, scratch_dir){
    res = foreach(i = 1:length(myrxcui), .combine = rbind )%dopar%{
        if (i %% 100 == 0 ) cat("\r",i," of ", length(myrxcui),"\t\t")
        term = myrxcui[i]
        url = paste("https://rxnav.nlm.nih.gov/REST", paste('/approximateTerm?term=',term,'&maxEntries=1&option=1',sep=""), sep="")
        url = gsub(" ", "%", url)
        destfile = paste(scratch_dir,"/approximateTerm",i,".xml",sep="")
        download.file(url, destfile = destfile, quiet=TRUE)
        check = try (x <- xmlParse(read_xml(destfile)), silent=T)
        if (!inherits(check, "try-error")){
            ns = getNodeSet(x,'/rxnormdata/approximateGroup/candidate', fun = xmlToList)
            ns = as.data.table(do.call("rbind", ns))
            if(nrow(ns)>0){
                ns[, input:=term]
                ns 
            }
        }
    }
    res[, rxcui := as.integer(rxcui) ]
    res[, rxaui := as.integer(rxaui) ]
    res[, score := as.integer(score) ]
    res[, input := as.character(input) ]
    res[, rank := NULL ]
    return(res)
}

eval_partial_matches = function(pmatches){
    set.seed(1984)
    torvw = c()
    for (x in unique(pmatches$algorithm)){
        i = 0.5
        j = 0.6
        bin_name = 0
        while (j < 1) {
            tmp = pmatches[alreadymapped==0 & algorithm==x & (score>i & score<=j)]
            print(nrow(tmp))
            if (nrow(tmp)>0){
                if (nrow(tmp)>10){
                    tmp = tmp[sample(nrow(tmp), 10)]
                } else {
                    tmp = tmp[sample(nrow(tmp), nrow(tmp))]
                }
                tmp[,bin:=bin_name]   
                torvw = rbind(torvw, tmp)
            }
            i = j
            j = j + 0.1
            bin_name = bin_name + 1
        }
    }
    torvw = torvw[sample(nrow(torvw),nrow(torvw))]
    
    ##manually review these entries to generate a confidence level in different match score bins
    for (row in 1:nrow(torvw)){
        cat(row,' out of ',nrow(torvw),'\n')
        cat(unlist(torvw[row,2]),'\n',unlist(torvw[row,3]),'\n')
        decision = as.integer(as.character(readLines(con='stdin',1)))
        torvw [ row, approve := decision ]
        torvw [ row+1, approve := "LEFT_OFF_HERE" ]
        save(torvw, file="tmp.Rdata")
        cat('\n')
    }
    torvw[ is.na(approve), approve := 0 ] 
    save(torvw, file="tmp_manual_review.Rdata") 
      
    ##format matrix for plotting
    mtx = matrix(NA, nrow=25, ncol=7)
    row = 1
    for (i in unique(torvw$algorithm)){
        thresh = 0.5
        for (j in 1:5) {
            mtx[row, 1] = i 
            mtx[row, 2] = thresh
            mtx[row, 3] = nrow(pmatches[algorithm==i & alreadymapped==0 & score >= thresh])
            mtx[row, 4] = nrow(torvw[algorithm==i  & alreadymapped==0 & score >= thresh])
            mtx[row, 5] = nrow(torvw[algorithm==i  & alreadymapped==0 & score >= thresh & approve==1])
            mtx[row, 6] = nrow(torvw[algorithm==i  & alreadymapped==0 & score >= thresh & approve==1])/nrow(torvw[algorithm==i & alreadymapped==0 &score >= thresh])
            mtx[row, 7] = nrow(pmatches[algorithm==i & alreadymapped==0 & score >= thresh])/nrow(pmatches[algorithm==i & alreadymapped==0])
            row = row + 1 
            thresh = thresh + 0.1
        }
    }
    mtx = as.data.table(mtx)[,.(algorithm=V1,threshold=as.numeric(V2),
                                total=as.numeric(V3),n_checked=as.numeric(V4),n_correct=as.numeric(V5),
                                pct_correct=as.numeric(V6),pct_remain=as.numeric(V7), pct_lost=1-as.numeric(V7))]
    pdf("~/www/figures/rx/partial_string_matching_performance_summary.pdf", height=8.5, width=11)
    ggplot(mtx,aes(pct_correct,pct_remain)) +
        geom_point(aes(col=algorithm), size=3) +
        geom_line(aes(col = algorithm), linetype = 2) +
        geom_text_repel(data = subset(mtx, threshold>0.5),aes(label=threshold, col = algorithm), fontface = 'bold', segment.color = "grey50") + 
        xlim(0,1) + ylim(0,1) + xlab("Percent Matches Correct") + ylab("Percent Matches Remain") + 
        ggtitle("Performance of 5 partial string matching algorithms")
    dev.off()
}

get_rxcui_family = function(my_status, scratchdir){
    active = unique(my_status[!is.na(remapped_rxcui)]$remapped_rxcui)
    family = foreach ( i = 1:length(active),.combine=rbind)%dopar%{
        if (i %% 100 == 0 ) cat("\r",i," of ", length(active),"\t\t")
        url = paste("https://rxnav.nlm.nih.gov/REST/rxcui/",active[i],"/status",sep="")
        destfile = paste(scratchdir, "/",active[i],".xml",sep="")
        download.file(url, destfile = destfile, quiet=TRUE)
        x = xmlParse(read_xml(destfile))
        ns = getNodeSet(x,'/rxnormdata/allRelatedGroup/conceptGroup/conceptProperties', fun = xmlToList)
        c1 = active[i]
        m = c()
        if (length(ns) > 0){
            for ( j in 1:length(ns)) {
                c2 = ns[[j]]$rxcui
                c3 = ns[[j]]$name
                c4 = ns[[j]]$umlscui
                c5 = ns[[j]]$tty
                m = rbind(m, as.data.table(matrix(c(c1,c2,c3,c4,c5),nrow=1,ncol=5)))
            }
            colnames(m) = c("remapped_rxcui","rel_rxcui","rel_name","rel_umlscui","rel_tty")
        }
        m
    }
    family = merge(family, family[remapped_rxcui==rel_rxcui,.(remapped_rxcui,remapped_name=rel_name,remapped_umlscui=rel_umlscui,remapped_tty=rel_tty)],
                   by="remapped_rxcui")
    family = family[,.(rxcui1=remapped_rxcui, rxcui2=rel_rxcui,name1=remapped_name,name2=rel_name,tty1=remapped_tty,tty2=rel_tty)]
    return(family)
}

rmsalts = function(names){
    words = unique(unlist(strsplit(names, split=" ")))
    salts = words[grep("ate$", words)]
    salts = c(salts,words[ words %in% c("sodium", "chloride", "bromide", "iodide","potassium", "choride", "magnesium", "calcium", "hydrabamine", "hydrabramine", "fosfatex", "zinc") ])
    new.salts = c()
    new.words = words[!words %in% salts]
    for(salt in salts) new.salts = c(new.salts, new.words[grep(salt,new.words)])
    salts = c(salts, new.salts)
    salts = c(salts,c("extract", "fruit", "seed", "root", "bark", "leaf", "herb", "oil", "essential oil",
                      "spp", "spore", "powder", "acid", "syrup", "salt", "salts", "vaccine", "dry"))
    lst = strsplit(names, split=" ")
    for (salt in salts) lst = lapply( lst, function(x) x = x[x!=salt] )
    names = unlist( lapply(lst, paste, collapse=" ") )
    return(names)
}

get_rxcui_in = function (rxcui_list, scratchdir, use.foreach=TRUE){
    if (use.foreach){
        res = foreach ( i = 1:length(rxcui_list),.combine=rbind)%dopar%{
            if (i %% 100 == 0 ) cat("\r",i," of ", length(rxcui_list),"\t\t")
            url = paste("https://rxnav.nlm.nih.gov/REST/rxcui/",rxcui_list[i],"/related?tty=IN",sep="")
            dest = paste(scratchdir, "/", rxcui_list[i], ".xml", sep="")
            download.file(url, destfile = dest, quiet=TRUE)
            x = xmlParse(read_xml(dest))
            file.remove(dest)
            out = getNodeSet(x,'/rxnormdata/relatedGroup/conceptGroup/conceptProperties', fun = xmlToList)
            if (length(out) > 0) {
                out = paste(unlist(lapply( out, function (x) x$rxcui)), collapse="|")
                out = data.table( rxcui_input = rxcui_list[i], rxcui_output = out )
                out
            }
        }
        return(res)
    } else {
        res = c()
        for ( i in 1:length(rxcui_list)) {
            cat("\r",i," of ", length(rxcui_list),"\t\t")
            url = paste("https://rxnav.nlm.nih.gov/REST/rxcui/",rxcui_list[i],"/related?tty=IN",sep="")
            dest = paste(scratchdir, "/", rxcui_list[i], ".xml", sep="")
            download.file(url, destfile = dest, quiet=TRUE)
            x = xmlParse(read_xml(dest))
            file.remove(dest)
            out = getNodeSet(x,'/rxnormdata/relatedGroup/conceptGroup/conceptProperties', fun = xmlToList)
            if (length(out) > 0) {
                out = paste(unlist(lapply( out, function (x) x$rxcui)), collapse="|")
                out = data.table( rxcui_input = rxcui_list[i], rxcui_output = out )
                res = rbind(res, out)
            }
        }
        return(res)               
    }
}

get_rxcui_status = function (rxcui_list, scratchdir){
    res <- c()
    ##res = foreach ( i = 1:length(rxcui_list),.combine=rbind)%dopar%{
    for ( i in 1:length(rxcui_list)){
        if (i %% 100 == 0 ) cat("\r",i," of ", length(rxcui_list),"\t\t")
        url = paste("https://rxnav.nlm.nih.gov/REST/rxcui/",rxcui_list[i],"/status",sep="")
        destfile = paste(scratchdir, "/", rxcui_list[i],".xml", sep="")
        download.file(url, destfile = destfile, quiet=TRUE)
        x = xmlParse(read_xml(destfile))
        file.remove(destfile)
        ns = getNodeSet(x,'/rxnormdata', fun = xmlToList)
        c1 = rxcui_list[i]
        c2 = ns[[1]]$rxcuiStatus$status
        c3 = c4 = c5 = c6 = NA
        if ( length(ns[[1]]$rxcuiStatus) > 2 ) {
            c3 = ns[[1]]$rxcuiStatus$minConceptGroup$minConcept$rxcui
            c4 = ns[[1]]$rxcuiStatus$minConceptGroup$minConcept$name
            c5 = ns[[1]]$rxcuiStatus$minConceptGroup$minConcept$tty
            if (c2 == "Remapped" ) {
                c6 = ns[[1]]$rxcuiStatus$remappedDate}
        }
        m = as.data.table(matrix(c(c1,c2,c3,c4,c5,c6),nrow=1,ncol=6))
        colnames(m) = c("mapped_rxcui","mapped_rxcui_status","remapped_rxcui","remapped_name","remapped_tty","remapped_date")
        ##m
        res <- rbind(res,m)
    }
    return(res)
}

load_ddx = function(path, date) {

    ##arguments are path to directory containing release and date of release in lower case (eg "sep_1_2017")
    ##helper file needs to be made manually from the pdf describing file structure (should be consistent mostly between releases)

    setwd(path)

    ##files
    main = paste("main_files_export_", date, sep="")
    help = paste("subsidiary_files_export_", date, sep="")
    whod = paste("who_ddx_", date, sep="")

    if (file.exists(paste(whod,"/ddx.Rdata", sep=""))){
        cat("ddx.Rdata exists already, loading from directory", paste(path,"/",whod,sep=""), "\n")
        load(paste(whod,"/ddx.Rdata", sep=""))
        return(ddx)
    } else {
        cat("ddx.Rdata does not exist, creating\n")
        mainfiles = list.files(main, full.names=T)
        names(mainfiles) = gsub(".txt", "", tstrsplit(mainfiles, split="/", keep=2L)[[1]])
        helpfiles = list.files(help, full.names=T)
        names(helpfiles) = gsub(".txt", "", tstrsplit(helpfiles, split="/", keep=2L)[[1]])
        info = fread("helper_file.txt", header=T)
        mainfiles.ddx = list.files(whod, full.names=T)
        names(mainfiles.ddx) = gsub(".txt", "", tstrsplit(mainfiles.ddx, split="/", keep=2L)[[1]])
        info.ddx = fread(paste(whod,"/helper_file.txt", sep=''))

        ##make
        ddx = list("key"=info.ddx)
        for (h in c("main", "subsid")){
            ddx[[h]] = list()
            for (i in unique(info.ddx[filegroup==h]$file)){
                ddx[[h]][[i]] = fread(mainfiles.ddx[i], sep='\t', header=F)
                for (j in unique(info.ddx[filegroup==h & file==i]$column)){
                    cat(h, i, j, "\n")
                    str = info.ddx[filegroup==h & file==i & column==j]$start
                    stp = info.ddx[filegroup==h & file==i & column==j]$stop
                    ddx[[h]][[i]][ , (j) := trimws(substr(V1, start = str, stop = stp))]
                    ddx[[h]][[i]][get(j) == "", (j):= NA]
                }
                ddx[[h]][[i]]$V1 = NULL
                ddx[[h]][[i]] = ddx[[h]][[i]][, colSums(is.na(ddx[[h]][[i]])) != nrow(ddx[[h]][[i]]), with=F]
            }
        }
        ddx$main$MP[, Drug_name.mod:=makemod(Drug_name)]
        ddx$subsid$SUN[, Substance_text.mod:=makemod(Substance_text)]
        
        ## Merge ddx ids into single matrix -----------------------------------------
        mm = merge( ddx$main$ING[,.(Ingredient_Id,Substance_Id,Medicinalprod_Id,Pharmproduct_Id)], ddx$subsid$SUN[,.(Substance_Id,Substance_text.mod, CAS_number)], all=TRUE )
        mm = merge( mm, ddx$main$PP[,.(Pharmproduct_Id,Medicinalprod_Id,Number_of_ingredients)], by = c("Pharmproduct_Id","Medicinalprod_Id"), all=TRUE )
        mm = merge(mm, ddx$main$MP[,.(Medicinalprod_Id,Medicinalprod_text.mod=Drug_name.mod)], by="Medicinalprod_Id", all=TRUE)
        mm[, Number_of_ingredients:= as.numeric(Number_of_ingredients) ]
        x = as.data.table(table(mm$Medicinalprod_Id))[,.(Medicinalprod_Id=V1, N)]
        tmp = merge(x, unique(mm[,.(Medicinalprod_Id, Number_of_ingredients)]), by="Medicinalprod_Id")
        nrow(tmp[ N !=Number_of_ingredients ]) # 0 ... good, sanity check
        ddx$mystuff = list()
        ddx$mystuff$idmap = mm[,.(Medicinalprod_Id, Pharmproduct_Id, Substance_Id, Ingredient_Id, Medicinalprod_text.mod, Number_of_ingredients, Substance_text.mod, CAS_number)]
        
        ## Standardize drug names so drugs with same ingreds have same name ----------------------
        gt1 =  unique(ddx$mystuff$idmap[Number_of_ingredients>1]$Medicinalprod_Id) #all mp with >1 sub
        gt1.dt = ddx$mystuff$idmap[Number_of_ingredients>1]
        gt1.dt = merge(gt1.dt, gt1.dt[, list( indices = paste(unique(sort(as.numeric(Substance_Id))), collapse="_") ), by = list( Medicinalprod_Id) ], by="Medicinalprod_Id")
        lst = list()
        for (h in 1:length(gt1)){
            i = gt1[h]
            cat('\r', h, '\t\t')
            indices = paste(unique(sort(as.numeric(ddx$mystuff$idmap[Medicinalprod_Id==i]$Substance_Id))), collapse="_")
            if (indices %in% names(lst)){
                lst[[indices]] = c(lst[[indices]], i)
            } else {
                lst[[indices]] = i
            }
        }
        tmp = data.table( nmp=unlist(lapply(lst, length)), indices=names(unlist(lapply(lst, length))))
        gt1.dt = merge( gt1.dt, tmp, by="indices")
        gt1.dt = gt1.dt[nmp>1]
        gt1.dt = merge(gt1.dt, gt1.dt[, list(mp = paste(unique(sort(as.integer(Medicinalprod_Id))), collapse="|"),
                                             mpt = paste(unique(sort(Medicinalprod_text.mod)), collapse="|") ), by = list(indices) ], by="indices")
        gt1.dt[,Medicinalprod_Id.unify:=mp]
        gt1.dt[,Medicinalprod_text.mod.unify:=mpt]
        gt1.dt[,mp:=NULL]
        gt1.dt[,mpt:=NULL]
        ddx$mystuff$idmap = merge(ddx$mystuff$idmap, gt1.dt[,.(Medicinalprod_Id,Substance_Id,Medicinalprod_Id.unify,Medicinalprod_text.mod.unify)], by=c("Medicinalprod_Id","Substance_> "), all=TRUE)
        ddx$mystuff$idmap[is.na(Medicinalprod_Id.unify), Medicinalprod_Id.unify:=Medicinalprod_Id]
        ddx$mystuff$idmap[Medicinalprod_Id.unify==Medicinalprod_Id, Medicinalprod_text.mod.unify:=Medicinalprod_text.mod]
        problems = as.data.table(table(ddx$mystuff$idmap$Medicinalprod_text.mod.unify))[N>1]$V1
        fixthese = unique(ddx$mystuff$idmap[ Medicinalprod_text.mod.unify %in% problems & Number_of_ingredients==1]$Medicinalprod_text.mod.unify) #single ingredient drugs with multiple med prod ids
        cnt = 1
        for (i in fixthese){
            cat("\r",cnt/length(fixthese),"\t\t")
            pick = paste(ddx$mystuff$idmap[Medicinalprod_text.mod.unify==i]$Medicinalprod_Id.unify, collapse="|")
            ddx$mystuff$idmap[Medicinalprod_text.mod.unify==i, Medicinalprod_Id.unify:=pick]
            cnt = cnt+1
        }
        
        tmp = unique(ddx$mystuff$idmap[Number_of_ingredients==1,.(Substance_text.mod,Medicinalprod_Id,Medicinalprod_text.mod,
                                                                  Number_of_ingredients)])
        tmp = tmp[, list(mp = paste(unique(sort(as.integer(Medicinalprod_Id))), collapse="|"),
                         mpt = paste(unique(sort(Medicinalprod_text.mod)), collapse="|")), by = list(Substance_text.mod,Number_of_ingredients) ]
        ddx$mystuff$idmap = merge(ddx$mystuff$idmap, tmp, by=c("Substance_text.mod","Number_of_ingredients"), all.x=TRUE)
        ddx$mystuff$idmap[!is.na(mp), Medicinalprod_Id.unify:=mp]
        ddx$mystuff$idmap[!is.na(mpt), Medicinalprod_text.mod.unify:=mpt]
        ddx$mystuff$idmap[, mp:=NULL]
        ddx$mystuff$idmap[, mpt:=NULL]
        tmp = unique(ddx$mystuff$idmap[Number_of_ingredients>1,.(Substance_text.mod,Substance_Id,
                                                                 Medicinalprod_Id.unify, Number_of_ingredients)])
        tmp = tmp[, list(Substance_text.mod.list = paste(unique(sort(Substance_text.mod)), collapse="|"),
                         Substance_Id.list = paste(unique(sort(as.integer(Substance_Id))), collapse="|")),
                  by = list(Medicinalprod_Id.unify,Number_of_ingredients) ]
        ddx$mystuff$idmap <-  merge(ddx$mystuff$idmap, tmp, by=c("Medicinalprod_Id.unify","Number_of_ingredients"), all.x=TRUE)
        ddx$mystuff$idmap[Number_of_ingredients==1, Substance_text.mod.list:=Substance_text.mod]
        ddx$mystuff$idmap[Number_of_ingredients==1, Substance_Id.list:=Substance_Id]
        uniqueN(ddx$mystuff$idmap[,.(Medicinalprod_text.mod.unify,Medicinalprod_Id.unify)]) #[1] 38716
        uniqueN(ddx$mystuff$idmap[,.(Medicinalprod_text.mod.unify)]) #[1] 38716
        uniqueN(ddx$mystuff$idmap[,.(Medicinalprod_Id.unify)]) #[1] 38716
        uniqueN(ddx$mystuff$idmap[,.(Substance_text.mod.list,Substance_Id.list)]) #[1] 37826
        uniqueN(ddx$mystuff$idmap[,.(Substance_text.mod.list)]) #[1] 37826
        uniqueN(ddx$mystuff$idmap[,.(Substance_Id.list)]) #[1] 37826
        mytmp = data.table(Substance_text.mod=unique(ddx$mystuff$idmap$Substance_text.mod),
                           Substance_text.mod.trim=rmsalts(unique(ddx$mystuff$idmap$Substance_text.mod)))
        mytmp[Substance_text.mod.trim=="",Substance_text.mod.trim:=NA]
        ddx$mystuff$idmap = merge(ddx$mystuff$idmap, mytmp, by="Substance_text.mod")
        save(ddx,file=paste(whod,"/ddx.Rdata", sep=""))
        return(ddx)
    }
}
        

load_vigibase = function(path, date, filelist) {

    cat("Retrieving VigiBase data files from", path, "\n")
    real_start = Sys.time()
    cat("Start time is", as.character(real_start), "\n")    
    setwd(path)

    ## directories
    main = paste("main_files_export_", date, sep="")
    help = paste("subsidiary_files_export_", date, sep="")
    whod = paste("who_ddx_", date, sep="")
    
    ## files
    mainfiles = list.files(main, full.names=T)
    names(mainfiles) = gsub(".txt", "", tstrsplit(mainfiles, split="/", keep=2L)[[1]])
    helpfiles = list.files(help, full.names=T)
    names(helpfiles) = gsub(".txt", "", tstrsplit(helpfiles, split="/", keep=2L)[[1]])
    info = fread("helper_file.txt", header=T)
    mainfiles.ddx = list.files(whod, full.names=T)
    names(mainfiles.ddx) = gsub(".txt", "", tstrsplit(mainfiles.ddx, split="/", keep=2L)[[1]])
    info.ddx = fread(paste(whod,"/helper_file.txt", sep=''))

    ## initiate output
    ret.list = list()
    
    ## demo
    if ("demo" %in% filelist){
        start = Sys.time()
        cat("LOADING DEMO\n")
        demo = fread(mainfiles["DEMO"], sep="\t", header=F)
        demo = main2subsid(info=info, dt=demo, dtname="DEMO", helpfiles=helpfiles)
        ret.list[["demo"]] = demo                     
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }
    
    ## adr
    if ("adr" %in% filelist){
        start = Sys.time()
        cat("LOADING ADR\n")
        adr = fread(mainfiles["ADR"], sep="\t", header=F)
        adr = main2subsid(info=info, dt=adr, dtname="ADR", helpfiles=helpfiles)
        ret.list[["adr"]] = adr
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }
    
    ## drug
    if ("drug" %in% filelist){
        start = Sys.time()
        cat("LOADING DRUG\n")
        drug = fread(mainfiles["DRUG"], sep="\t", header=F)
        drug = main2subsid(info=info, dt=drug, dtname="DRUG", helpfiles=helpfiles)
        ret.list[["drug"]] = drug
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }
    
    ## ind
    if ("ind" %in% filelist){
        start = Sys.time()
        cat("LOADING IND\n")
        ind = fread(mainfiles["IND"], sep="\t", header=F)
        ind = main2subsid(info=info, dt=ind, dtname="IND", helpfiles=helpfiles)
        ret.list[["ind"]] = ind
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }
    
    ## out
    if ("out" %in% filelist){
        start = Sys.time()
        cat("LOADING OUT\n")
        out = fread(mainfiles["OUT"], sep="\t", header=F)
        out = main2subsid(info=info, dt=out, dtname="OUT", helpfiles=helpfiles)
        ret.list[["out"]] = out
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }
    
    ## link
    if ("link" %in% filelist){
        start = Sys.time()
        cat("LOADING LINK\n")
        link = fread(mainfiles["LINK"], sep="\t", header=F)
        link = main2subsid(info=info, dt=link, dtname="LINK", helpfiles=helpfiles)
        ret.list[["link"]] = link
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }
    
    ## srce
    if ("srce" %in% filelist){
        start = Sys.time()
        cat("LOADING SRCE\n")
        srce = fread(mainfiles["SRCE"], sep="\t", header=F)
        srce = main2subsid(info=info, dt=srce, dtname="SRCE", helpfiles=helpfiles)
        ret.list[["srce"]] = srce
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }

    ## followup
    if ("followup" %in% filelist){
        start = Sys.time()
        cat("LOADING FOLLOWUP\n")
        followup = fread(mainfiles["FOLLOWUP"], sep="\t", header=F)
        followup = main2subsid(info=info, dt=followup, dtname="FOLLOWUP", helpfiles=helpfiles)
        ret.list[["followup"]] = followup
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }
        
    ## drug dictionary (ddx)
    if ("ddx" %in% filelist){
        start = Sys.time()
        cat("LOADING DDX\n")
        ddx = load_ddx(path, date)
        ret.list[["ddx"]] = ddx
        end = Sys.time()
        cat(paste("Time elapsed:",  round(as.numeric(difftime(end, start, units="mins")), 2), "minutes"), "\n")
        cat("\n")
    }

    real_end = Sys.time()
    cat("End time is", as.character(real_end), "\n")    
    cat(paste("Total time elapsed:",  round(as.numeric(difftime(real_end, real_start, units="mins")), 2), "minutes"), "\n")
    cat("DONE\n")
    return(ret.list)

}
 

fixmupt = function(myname, ptlist, mdr.wide, start="pt"){
    if (start == "pt"){
        myhlt = unique(mdr.wide[primary_soc_fg=="Y" & pt_code %in% ptlist]$hlt_code)
        myhlgt = unique(mdr.wide[primary_soc_fg=="Y" & pt_code %in% ptlist]$hlgt_code)
        mysoc = unique(mdr.wide[primary_soc_fg=="Y" & pt_code %in% ptlist]$soc_code)
        if(length(myhlt)==1){
            mydesc = unique(mdr.wide[hlt_code==myhlt]$hlt_name)
            out = data.table( inputcode = myname, lev = "hlt", meddra = myhlt, meddra_name=mydesc)
        } else if (length(myhlgt)==1) {
            mydesc = unique(mdr.wide[hlgt_code==myhlgt]$hlgt_name)
            out = data.table( inputcode = myname, lev = "hlgt", meddra = myhlgt, meddra_name=mydesc)
        } else if (length(mysoc)==1) {
            mydesc = unique(mdr.wide[soc_code==mysoc]$soc_name)
            out = data.table( inputcode = myname, lev = "soc", meddra = mysoc, meddra_name=mydesc)
        } else {
            mychoice = sample(ptlist, 1)
            mydesc = unique(mdr.wide[pt_code==mychoice]$pt_name)
            out = data.table( inputcode = myname, lev = "pt", meddra = mychoice, meddra_name=mydesc)
        }
    }
    if (start == "hlt"){
        myhlgt = unique(mdr.wide[primary_soc_fg=="Y" & hlt_code %in% ptlist]$hlgt_code)
        mysoc = unique(mdr.wide[primary_soc_fg=="Y" & hlt_code %in% ptlist]$soc_code)
        if (length(myhlgt)==1) {
            mydesc = unique(mdr.wide[hlgt_code==myhlgt]$hlgt_name)
            out = data.table( inputcode = myname, lev = "hlgt", meddra = myhlgt, meddra_name=mydesc)
        } else if (length(mysoc)==1) {
            mydesc = unique(mdr.wide[soc_code==mysoc]$soc_name)
            out = data.table( inputcode = myname, lev = "soc", meddra = mysoc, meddra_name=mydesc)
        } else {
            mychoice = sample(ptlist, 1)
            mydesc = unique(mdr.wide[hlt_code==mychoice]$hlt_name)
            out = data.table( inputcode = myname, lev = "hlt", meddra = mychoice, meddra_name=mydesc )
        }
    }
    if (start == "hlgt"){
        mysoc = unique(mdr.wide[primary_soc_fg=="Y" & hlt_code %in% ptlist]$soc_code)
        if (length(mysoc)==1) {
            mydesc = unique(mdr.wide[soc_code==mysoc]$soc_name)
            out = data.table( inputcode = myname, lev = "soc", meddra = mysoc, meddra_name=mydesc)
        } else {
            mychoice = sample(ptlist, 1)
            mydesc = unique(mdr.wide[hlgt_code==mychoice]$hlgt_name)
            out = data.table( inputcode = myname, lev = "hlgt", meddra = mychoice, meddra_name=mydesc )
        }
    }
    return(out)
}

disprostats = function(dt) {

    if ( !identical( sort(colnames(dt)), c("ae","drug","report") ) ) stop("ERROR: colnames for input must be 'report' 'drug' and 'ae'")        
    drugs = unique(dt$drug)
    aes = unique(dt$ae)
    
    ##contingency table
    cat("Making contingency table...")
    nreports = uniqueN(dt[,.(report)])
    dt = unique(dt[,.(report,drug,ae)])    
    tmp0 = dt[, .N, by=list(drug,ae)][,.(drug,ae,a=N)]
    tmp1 = CJ(drugs, aes, unique = TRUE)[,.(drug=drugs, ae=aes)]
    tmp2 = dt[, .N, by=list(drug)][,.(drug,n.withDRUG=N)]
    tmp3 = dt[, .N, by=list(ae)][,.(ae,n.withAE=N)]
    contab = merge(tmp0, tmp1, by=c("drug","ae"), all.y=TRUE)
    contab[is.na(a), a:=0]
    contab = merge(contab, tmp2, by="drug")
    contab = merge(contab, tmp3, by="ae")
    contab[, b := (n.withDRUG - a) ]
    contab[, c := (n.withAE - a) ]
    contab[, d := nreports - (a + b + c) ]
    contab = contab[,.(drug, ae, a, b, c, d)]
    cat("done\n")

    ##ROR
    cat("Calculating ROR stats...")
    contab[ , ror := (a/c)/(b/d) ]
    contab[ ror > 0, ror_se := sqrt( (1/a) + (1/b) + (1/c) + (1/d) ) ]
    contab[ ror > 0, ror_95ci_ul := exp(1)^(log(ror) + (1.96*ror_se)) ]
    contab[ ror > 0, ror_95ci_ll := exp(1)^(log(ror) - (1.96*ror_se)) ]
    cat("done\n")
    
    ##PRR
    cat("Calculating PRR stats...")
    contab[, prr := (a/(a+b))/(c/(c+d)) ]
    contab[ prr > 0, prr_se := sqrt( (1/a) - (1/(a+b)) + (1/c) + (1/(c+d)) ) ]
    contab[ prr > 0, prr_95ci_ul := exp(1)^(log(prr) + (1.96*prr_se)) ]
    contab[ prr > 0, prr_95ci_ll := exp(1)^(log(prr) - (1.96*prr_se)) ]
    cat("done\n")
        
    ##INFORMATION COMPONENT
    cat("Calculating IC stats...")
    tmp4 = copy(contab)
    tmp4[, C := a+b+c+d ]   #total reports in database
    tmp4[, c_i := a+b] #nreports of drug i in database
    tmp4[, c_j := a+c] #nreports of ae j in database
    tmp4[, c_ij := a ] #nreports of drug and ae in database
    gamma_ij   = 1
    alpha_i    = 1
    alpha      = 2
    beta_j     = 1
    beta       = 2
    tmp4[, gamma := gamma_ij * ((C + alpha)/(c_i+alpha_i)) * ((C + beta)/(c_j+beta_j)) ]
    tmp4[, IC := log2(((c_ij + gamma_ij)*(C+alpha)*(C+beta)) / ((C + gamma)*(c_i+alpha_i)*(c_j+beta_j))) ]
    tmp4[, IC_var := (((C-c_ij+gamma-gamma_ij)/((c_ij+gamma_ij)*(1+C+gamma)))+((C-c_i+alpha-alpha_i)/((c_i+alpha_i)*(1+C+alpha)))+((C-c_j+beta-beta_j)/((c_j+beta_j)*(1+C+beta))))/log(2)^2 ]
    tmp4[, IC_sd := sqrt(IC_var) ]
    contab = merge(contab, tmp4[,.(drug,ae,IC,IC_var,IC_sd)], by=c("drug","ae")) 
    cat("done\n")
    
    ##CHI SQ with Yates correction
    cat("Calculating Chi-Squared (with Yates' correction) stats...")
    mychistat = function (x) chisq.test( matrix(x, ncol=2, byrow=TRUE), correct=T )$statistic
    mychip = function (x) chisq.test( matrix(x, ncol=2, byrow=TRUE), correct=T )$p.value
    tmp5 = contab[a>0, list( chisq.yates = mychistat( c(a,b,c,d) ) ), by=list(drug,ae) ]
    tmp6 = contab[a>0, list( chisq.yates = mychip( c(a,b,c,d) ) ), by=list(drug,ae) ]
    tmp7 = merge(tmp5, tmp6, by=c("drug","ae"), suffixes=c(".stat",".p"))    
    contab = merge(contab, tmp7, all.x=TRUE)
    cat("done\n")
    
    ##rename cols
    colnames(contab)[3] = "nreport.with.drug.with.ae"
    colnames(contab)[4] = "nreport.with.drug.without.ae"
    colnames(contab)[5] = "nreport.with.ae.without.drug"
    colnames(contab)[6] = "nreport.without.drug.without.ae"
        
    ##return
    return(contab)
    
}

disprostats.check <- function(input.data, steps=1000) {
    ## Using cutoffs for prr/ror from Table 1 of Puijenbroek et al 2003
    ##
    ## NOTE: these cutoffs were not strictly used to define signal vs no signal, but rather to
    ##       assess specificity, sensitivity, PPV and NPV of the PRR and ROR (using a different test as
    ##       as the standard to define ground truth)
    ncombos = uniqueN(input.data[,.(drug,ae)])
    ncombos.nosdr = uniqueN(input.data[sider==0,.(drug,ae)])
    ncombos.sdr = uniqueN(input.data[sider==1,.(drug,ae)])
    check = c()
    for ( i in 1:steps) {
        print(i)
        rordat = input.data[ (ror-(1.96*ror_se))>1 & nreport.with.drug.with.ae>=i,]
        prrdat = input.data[ (prr-(1.96*prr_se))>1 & nreport.with.drug.with.ae>=i,]
        chidat = input.data[ chisq.yates.p.adj.bon<0.05 & chisq.yates.stat>4 & nreport.with.drug.with.ae>=i,]
        iccdat  = input.data[ (IC - (2*IC_sd))>0 & nreport.with.drug.with.ae>=i,]
        rortab = table(rordat$sider)
        prrtab = table(prrdat$sider)
        chitab = table(chidat$sider)
        icctab = table(icdat$sider)
        rorcnt = nrow(rordat)
        prrcnt = nrow(prrdat)
        chicnt = nrow(chidat)
        icccnt = nrow(icdat)
        rorcnt.nosdr = rortab[1]
        prrcnt.nosdr = prrtab[1]
        chicnt.nosdr = chitab[1]
        icccnt.nosdr = ictab[1]
        rorcnt.sdr = rortab[2]
        prrcnt.sdr = prrtab[2]
        chicnt.sdr = chitab[2]
        icccnt.sdr = ictab[2]
        rorfrac = rorcnt/ncombos
        prrfrac = prrcnt/ncombos
        chifrac = chicnt/ncombos
        iccfrac = icccnt/ncombos
        rorfrac.nosdr = rorcnt.nosdr/ncombos.nosdr
        prrfrac.nosdr = prrcnt.nosdr/ncombos.nosdr
        chifrac.nosdr = chicnt.nosdr/ncombos.nosdr
        iccfrac.nosdr = icccnt.nosdr/ncombos.nosdr
        rorfrac.sdr = rorcnt.sdr/ncombos.sdr
        prrfrac.sdr = prrcnt.sdr/ncombos.sdr
        chifrac.sdr = chicnt.sdr/ncombos.sdr
        iccfrac.sdr  = icccnt.sdr/ncombos.sdr
        check = rbind(check,
                      data.table(min.reports=i, sider="all", stat="ror", nkeep=rorcnt, frackeep = rorfrac),
                      data.table(min.reports=i, sider="all", stat="prr", nkeep=prrcnt, frackeep = prrfrac),
                      data.table(min.reports=i, sider="all", stat="chi", nkeep=chicnt, frackeep = chifrac),
                      data.table(min.reports=i, sider="all", stat="icc", nkeep=icccnt, frackeep = iccfrac),
                      data.table(min.reports=i, sider="no.sdr", stat="ror", nkeep=rorcnt.nosdr, frackeep = rorfrac.nosdr),
                      data.table(min.reports=i, sider="no.sdr", stat="prr", nkeep=prrcnt.nosdr, frackeep = prrfrac.nosdr),
                      data.table(min.reports=i, sider="no.sdr", stat="chi", nkeep=chicnt.nosdr, frackeep = chifrac.nosdr),
                      data.table(min.reports=i, sider="no.sdr", stat="icc", nkeep=icccnt.nosdr, frackeep = iccfrac.nosdr),
                      data.table(min.reports=i, sider="sdr", stat="ror", nkeep=rorcnt.sdr, frackeep = rorfrac.sdr),
                      data.table(min.reports=i, sider="sdr", stat="prr", nkeep=prrcnt.sdr, frackeep = prrfrac.sdr),
                      data.table(min.reports=i, sider="sdr", stat="chi", nkeep=chicnt.sdr, frackeep = chifrac.sdr),
                      data.table(min.reports=i, sider="sdr", stat="icc", nkeep=icccnt.sdr, frackeep = iccfrac.sdr))
    }
    return(check)
}

mydat_stats = function ( mydat, mysdr, mymdi, mymdr ) {

    cat("vigibase drug/ae pairs\n")
    vgb.drugae = unique(mydat[, .(UMCReportId,Drug_Id,adr.vbid,meddra=adr.meddra,rxcui=rxcui_output)])

    cat("vigibase drug/ind pairs\n")
    vgb.drugind = unique(mydat[, .(UMCReportId,Drug_Id,meddra=ind.meddra,rxcui=rxcui_output)])

    cat("sider drug/ae pairs\n")
    sdr.drugae = unique(mysdr[meddra_types=="ae",.(rxcui=as.character(rxcui),meddra)]) #drug/ae pairs in sider

    cat("sider drug/ind pairs\n")
    sdr.drugind = unique(mysdr[meddra_types %in% c("ae|ind","ind"),.(rxcui=as.character(rxcui),meddra)]) #drug/ind pairs in sider

    cat("medi drug/ind pairs\n")
    mdi.drugind = unique(mymdi[,.(rxcui=as.character(rxcui),meddra)]) #drug/ind pairs in medi
    
    cat("all reports\n")
    all = unique(mydat[,.(UMCReportId)])

    cat("number of ae per report (for reports with >=1 ae)\n")
    tmp1 = unique(mydat[!is.na(adr.meddra.n),.(UMCReportId,n.ae=adr.meddra.n)])

    cat("number of drugs (and drug/ind pairs) per report\n")
    tmp2 = unique(mydat[,.(UMCReportId,Drug_Id,rxcui_output)])[,.N,by=UMCReportId][,.(UMCReportId,n.ing=N, n.ing.ind.pairs=N)]

    cat("number of drugs with no substance per report\n")
    tmp3 = unique(mydat[is.na(Substance_text.mod),.(UMCReportId,Drug_Id)])[,.N,by=UMCReportId][,.(UMCReportId,n.drug.with.no.sub=N)]

    cat("number of drugs mapped to rxnorm per report\n")
    tmp4 = unique(mydat[!is.na(rxcui_output),.(UMCReportId,Drug_Id,rxcui_output)])[,.N,by=UMCReportId][,.(UMCReportId,n.ing.mapped=N)]

    cat("number of  drugs with multiple ingredients per report\n")
    mymulti = unique(mydat[grep("|", rxcui_output, fixed=TRUE)]$rxcui_output) ##multi-ing drugs we found in mapping process
    tmp5a = unique(mydat[nSubstance>1,.(UMCReportId,Drug_Id)])[,.N,by=UMCReportId][,.(UMCReportId,n.drug.multi.ing.in.vb=N)]
    tmp5b = unique(mydat[rxcui_output %in% mymulti,.(UMCReportId,Drug_Id)])[,.N,by=UMCReportId][,.(UMCReportId,n.drug.multi.ing.by.map=N)]
    tmp5 = merge(tmp5a, tmp5b, all=TRUE)

    cat("number of indications mapped to meddra per report\n")
    tmp6 = unique(mydat[!is.na(ind.meddra),.(UMCReportId,ind.meddra)])[,.N,by=UMCReportId][,.(UMCReportId,n.ind.mapped=N)]

    cat("drug and ind mapped\n")
    tmp7 = unique(mydat[!is.na(ind.meddra)&!is.na(rxcui_output),.(UMCReportId, Drug_Id, rxcui_output)])[,.N,UMCReportId][,.(UMCReportId,n.ing.ind.pairs.2mapped=N)]

    cat("drug or ind mapped\n")
    tmp8 = unique(mydat[!is.na(ind.meddra)|!is.na(rxcui_output),.(UMCReportId, Drug_Id, rxcui_output)])[,.N,UMCReportId][,.(UMCReportId,n.ing.ind.pairs.1mapped=N)]

    cat("number of drugs suspected involved in ae per report\n")
    tmp9 = unique(mydat[Basis %in% c("Interacting", "Suspect"), .(UMCReportId,Drug_Id,rxcui_output)])[,.N,UMCReportId][,.(UMCReportId,n.ing.suspects=N)]

    cat("number of drug/ae pairs per report\n")
    tmp10 = vgb.drugae[,.N,UMCReportId][,.(UMCReportId,n.ing.ae.pairs=N)]

    cat("number of drug/ae pairs per report in sider\n")
    tmp11 = merge( sdr.drugae, vgb.drugae, by=c("rxcui","meddra") )[,.N,UMCReportId][,.(UMCReportId,n.ing.ae.pairs.in.sider=N)]

    cat("number of drug/ind pairs per report in sider\n")
    tmp12 = merge( sdr.drugind, vgb.drugind, by=c("rxcui","meddra") )[,.N,UMCReportId][,.(UMCReportId,n.ing.ind.pairs.in.sider=N)]

    cat("number of drug/ind pairs per report in medi\n")
    tmp13 = merge( mdi.drugind[!is.na(rxcui) & !is.na(meddra)], vgb.drugind, by=c("rxcui","meddra") )[,.N,UMCReportId][,.(UMCReportId,n.ing.ind.pairs.in.medi=N)]

    cat("number of drug/ind pairs per report sider and medi\n")
    x = rbind(mdi.drugind[!is.na(rxcui) & !is.na(meddra)], sdr.drugind)[,.N,list(rxcui,meddra)][N>1,.(rxcui,meddra)]
    tmp14 = merge(x, vgb.drugind, by=c("rxcui","meddra") )[,.N,UMCReportId][,.(UMCReportId,n.ing.ind.pairs.in.sider.and.medi=N)]

    cat("number of drug/ind pairs per report sider or medi\n")
    x = rbind(mdi.drugind[!is.na(rxcui) & !is.na(meddra)], sdr.drugind)[,.N,list(rxcui,meddra)][N>=1,.(rxcui,meddra)]
    tmp15 = merge(x, vgb.drugind, by=c("rxcui","meddra") )[,.N,UMCReportId][,.(UMCReportId,n.ing.ind.pairs.in.sider.or.medi=N)]

    cat("number of ae that are also an ind in same report\n")
    tmp16 = unique(mydat[adr.meddra == ind.meddra,.(UMCReportId,adr.vbid)])[,.N,UMCReportId][,.(UMCReportId,n.ae.equal.ind.in.report=N)]

    cat("number of ae that are an ind for drugs in report somewhere else in vigibase\n")
    vgb.drugae.mapped = unique(mydat[!is.na(rxcui_output) & !is.na(adr.meddra),.(rxcui=rxcui_output,meddra=adr.meddra)])
    vgb.drugind.mapped = unique(mydat[!is.na(rxcui_output) & !is.na(ind.meddra),.(rxcui=rxcui_output,meddra=ind.meddra)])
    same = merge( vgb.drugae.mapped, vgb.drugind.mapped )
    tmp17 =  unique(merge(mydat, same[,.(rxcui_output=rxcui, adr.meddra=meddra)],
                          by=c("rxcui_output","adr.meddra"))[,.(UMCReportId,adr.vbid)])[,.N,UMCReportId][,.(UMCReportId,n.ae.equal.ind.in.other.reports=N)]

    cat("number of ae that are an ind for drugs in report in sider\n")
    x = merge(vgb.drugae.mapped, sdr.drugind, by=c("rxcui","meddra"))
    x = merge(mydat, x[,.(rxcui_output=rxcui, adr.meddra=meddra)], by=c("rxcui_output","adr.meddra") )
    x = unique(x[,.(UMCReportId,Drug_Id,adr.vbid)])
    tmp18 = x[,.N,UMCReportId][,.(UMCReportId,n.ing.ae.pairs.as.ing.ind.pairs.in.sider=N)]

    cat("number of ae that are an ind for drugs in report in medi\n")
    x = merge(vgb.drugae.mapped, mdi.drugind[!is.na(rxcui) & !is.na(meddra)], by=c("rxcui","meddra"))
    x = merge(mydat, x[,.(rxcui_output=rxcui, adr.meddra=meddra)], by=c("rxcui_output","adr.meddra") )
    x = unique(x[,.(UMCReportId,Drug_Id,adr.vbid)])
    tmp19 = x[,.N,UMCReportId][,.(UMCReportId,n.ing.ae.pairs.as.ing.ind.pairs.in.medi=N)]

    cat("number of ae that are an ind for drugs in report in sider and medi\n")
    x = rbind(mdi.drugind[!is.na(rxcui) & !is.na(meddra)], sdr.drugind)[,.N,list(rxcui,meddra)][N>1,.(rxcui,meddra)]
    x = merge(vgb.drugae.mapped, x, by=c("rxcui","meddra"))
    x = merge(mydat, x[,.(rxcui_output=rxcui, adr.meddra=meddra)], by=c("rxcui_output","adr.meddra") )
    x = unique(x[,.(UMCReportId,Drug_Id,adr.vbid)])
    tmp20 = x[,.N,UMCReportId][,.(UMCReportId,n.ing.ae.pairs.as.ing.ind.pairs.in.sider.and.medi=N)]

    cat("number of ae that are an ind for drugs in report in sider or medi\n")
    x = rbind(mdi.drugind[!is.na(rxcui) & !is.na(meddra)], sdr.drugind)[,.N,list(rxcui,meddra)][N>=1,.(rxcui,meddra)]
    x = merge(vgb.drugae.mapped, x, by=c("rxcui","meddra"))
    x = merge(mydat, x[,.(rxcui_output=rxcui, adr.meddra=meddra)], by=c("rxcui_output","adr.meddra") )
    x = unique(x[,.(UMCReportId,Drug_Id,adr.vbid)])
    tmp21 = x[,.N,UMCReportId][,.(UMCReportId,n.ing.ae.pairs.as.ing.ind.pairs.in.sider.or.medi=N)]

    cat("lowest meddra level connecting ae and ind in a report\n")
    tmp22 = mydat[!is.na(adr.meddra) & !is.na(ind.meddra),.(UMCReportId,adr.meddra,adr.meddra.level,ind.meddra,ind.meddra.level)]
    tmp22.1 = copy(mymdr)
    tmp22.2 = copy(mymdr)
    colnames(tmp22.1) = paste( "adr.meddra.", colnames(tmp22.1), sep="" )
    colnames(tmp22.2) = paste( "ind.meddra.", colnames(tmp22.2), sep="" )
    tmp22.1 = unique(tmp22.1[adr.meddra.primary_soc_fg=="Y",.(adr.meddra.pt_code, adr.meddra.hlt_code, adr.meddra.hlgt_code, adr.meddra.soc_code)])
    tmp22.2a = unique(tmp22.2[ind.meddra.primary_soc_fg=="Y",.(ind.meddra.pt_code, ind.meddra.hlt_code, ind.meddra.hlgt_code, ind.meddra.soc_code)])
    tmp22.2b = unique(tmp22.2[,.(ind.meddra.hlt_code, ind.meddra.hlgt_code, ind.meddra.soc_code,ind.meddra.primary_soc_fg)])
    tmp22.2b.prim = tmp22.2b[ind.meddra.primary_soc_fg=="Y"]$ind.meddra.hlt_code
    tmp22.2b.x = unique(tmp22.2b[ ind.meddra.primary_soc_fg=="Y" & ind.meddra.hlt_code %in% tmp22.2b.prim,
                                 .(ind.meddra.hlt_code,ind.meddra.hlgt_code,ind.meddra.soc_code)])
    tmp22.2b.y = unique(tmp22.2b[!ind.meddra.hlt_code %in% tmp22.2b.x$ind.meddra.hlt_code,.(ind.meddra.hlt_code,ind.meddra.hlgt_code,ind.meddra.soc_code)])
    tmp22.2b = rbind(tmp22.2b.x, tmp22.2b.y)
    tmp22.2c = unique(tmp22.2[,.(ind.meddra.hlgt_code, ind.meddra.soc_code, ind.meddra.primary_soc_fg)])
    tmp22.2c.prim = tmp22.2c[ind.meddra.primary_soc_fg=="Y"]$ind.meddra.hlgt_code
    tmp22.2c.x = unique(tmp22.2c[ ind.meddra.primary_soc_fg=="Y" & ind.meddra.hlgt_code %in% tmp22.2c.prim,.(ind.meddra.hlgt_code,ind.meddra.soc_code)])
    tmp22.2c.y = unique(tmp22.2c[!ind.meddra.hlgt_code %in% tmp22.2c.x$ind.meddra.hlgt_code,.(ind.meddra.hlgt_code,ind.meddra.soc_code)])
    tmp22.2c = rbind(tmp22.2c.x, tmp22.2c.y)
    tmp22.3 = merge(tmp22, tmp22.1, by.x="adr.meddra", by.y="adr.meddra.pt_code")
    tmp22.3a = tmp22.3[ ind.meddra.level == "pt" ]
    tmp22.3b = tmp22.3[ ind.meddra.level == "hlt" ]
    tmp22.3c = tmp22.3[ ind.meddra.level == "hlgt" ]
    tmp22.3d = tmp22.3[ ind.meddra.level == "soc" ]
    tmp22.3aa = merge(tmp22.3a, tmp22.2a, by.x="ind.meddra", by.y="ind.meddra.pt_code")
    tmp22.3bb = merge(tmp22.3b, tmp22.2b, by.x="ind.meddra", by.y="ind.meddra.hlt_code")
    tmp22.3cc = merge(tmp22.3c, tmp22.2c, by.x="ind.meddra", by.y="ind.meddra.hlgt_code")
    tmp22.3dd = tmp22.3d
    tmp22.3aa[ , adr.meddra.pt_code := adr.meddra ]
    tmp22.3aa[ , ind.meddra.pt_code := ind.meddra ]
    tmp22.3bb[ , adr.meddra.pt_code := adr.meddra ]
    tmp22.3bb[ , ind.meddra.pt_code := NA ]
    tmp22.3bb[ , ind.meddra.hlt_code :=ind.meddra ]
    tmp22.3cc[ , adr.meddra.pt_code := adr.meddra ]
    tmp22.3cc[ , ind.meddra.pt_code := NA ]
    tmp22.3cc[ , ind.meddra.hlt_code := NA ]
    tmp22.3cc[ , ind.meddra.hlgt_code :=ind.meddra ]
    tmp22.3dd[ , adr.meddra.pt_code := adr.meddra ]
    tmp22.3dd[ , ind.meddra.pt_code := NA ]
    tmp22.3dd[ , ind.meddra.hlt_code := NA ]
    tmp22.3dd[ , ind.meddra.hlgt_code := NA ]
    tmp22.3dd[ , ind.meddra.soc_code := ind.meddra ]
    tmp22 = rbind(tmp22.3aa, tmp22.3bb, tmp22.3cc, tmp22.3dd)
    tmp22[!is.na(adr.meddra.pt_code) & !is.na(ind.meddra.pt_code) & adr.meddra.pt_code==ind.meddra.pt_code, meddra.level.where.adr.and.ind.equate := 0 ]
    tmp22[!is.na(adr.meddra.hlt_code) & !is.na(ind.meddra.hlt_code) & adr.meddra.hlt_code==ind.meddra.hlt_code & is.na(meddra.level.where.adr.and.ind.equate),
          meddra.level.where.adr.and.ind.equate := 1 ]
    tmp22[!is.na(adr.meddra.hlgt_code) & !is.na(ind.meddra.hlgt_code) & adr.meddra.hlgt_code==ind.meddra.hlgt_code & is.na(meddra.level.where.adr.and.ind.equate),
          meddra.level.where.adr.and.ind.equate := 2 ]
    tmp22[!is.na(adr.meddra.soc_code) & !is.na(ind.meddra.soc_code) & adr.meddra.soc_code==ind.meddra.soc_code & is.na(meddra.level.where.adr.and.ind.equate),
          meddra.level.where.adr.and.ind.equate := 3 ]
    tmp22[meddra.level.where.adr.and.ind.equate==0, tmpcol := "pt" ]
    tmp22[meddra.level.where.adr.and.ind.equate==1, tmpcol := "hlt" ]
    tmp22[meddra.level.where.adr.and.ind.equate==2, tmpcol := "hlgt" ]
    tmp22[meddra.level.where.adr.and.ind.equate==3, tmpcol := "soc" ]
    tmp22[is.na(meddra.level.where.adr.and.ind.equate), tmpcol := "none" ]
    tmp22[, tmpcol := factor(tmpcol, levels = c("pt", "hlt", "hlgt", "soc", "none"), ordered=TRUE) ]
    tmp22 = tmp22[,list(lowest.meddra.level.where.adr.and.ind.equate = min(tmpcol)),UMCReportId]
    tmp22[lowest.meddra.level.where.adr.and.ind.equate=="none",lowest.meddra.level.where.adr.and.ind.equate:=NA]
    
    cat("combine\n")
    mylst = list( tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp18, tmp19, tmp20, tmp21, tmp22)
    for ( i in 1:length(mylst)){all = merge( all, mylst[[i]], by="UMCReportId", all.x=TRUE )}

    cat("make NAs 0 where appropriate\n")
    myfun = function(x) {x[is.na(x)] = 0}
    for(i in colnames(all)[grep("n.", colnames(all))]) all[[i]][is.na(all[[i]])] = 0

    return(all)
    
}


