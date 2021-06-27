# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ ME
# 
# SCZ SNPs come from: https://www.medrxiv.org/content/10.1101/2020.09.12.20192922v1 (supplementary table 3)
# Gene location and name data comes from: https://grch37.ensembl.org/biomart/martview/0820cab9ba1c71dc01be7f18d4db0fab
#
# NOTE: These results are not included in the final analyses and as such are not included in the raw_data files previously provided
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Store user home directory and navigate to it
user_home_directory <- readline(prompt="Enter path where you have saved genetic_and_pharmacologic_overlap_in_schizophrenia-main:")
setwd(user_home_directory)

# Get the official home directory of the user
home_dir <- getwd()

# Navigate to the raw data zip file and unzip
setwd(home_dir)

library(data.table)
options(width=180)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load gene symbols to ID mapping
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene_names <- fread("/Users/Liam/Desktop/all_gene_names.tsv",sep="\t")
setnames(gene_names, "Approved symbol", "symbol")
setnames(gene_names, "Ensembl gene ID", "ensembl_id")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load gene locations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene_locations <- fread("/Users/Liam/Desktop/all_gene_locations.txt",sep="\t")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create working dataset with gene symbol and location
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene_info <- merge(gene_names,gene_locations,by="ensembl_id")
gene_info[,symbol := toupper(symbol)]
gene_info[,gene_location := gene_start_bp]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load and format PGC3 SZ SNPs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
snps <- fread("/Users/Liam/Desktop/all_scz_snps_final.txt", sep="\t")[,.(chromosome=as.character(chromosome),snp=top_index,p_value=top_P,snp_location=top_pos)]
snps[chromosome == "23", chromosome := "X"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up for looping to identify closest gene to each snp
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
list_of_snps <- unique(snps$snp)
count <- 1
res <- c()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify closest gene to each SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in list_of_snps) {
	print(paste0("Analyzing SNP #",count," out of ",length(list_of_snps)))

	curr_snp <- i
	curr_p <- unique(snps[snp==curr_snp]$p_value)
	curr_pos <- unique(snps[snp==curr_snp]$snp_location)
	curr_chromosome <- unique(snps[snp==curr_snp]$chromosome)

	gene_info_copy <- copy(gene_info)
	gene_info_copy <- gene_info_copy[chromosome==curr_chromosome,]
	gene_info_copy[,snp_location := curr_pos]
	gene_info_copy[,distance_to_snp := abs(snp_location-gene_location)]

	closest_gene_distance <- min(gene_info_copy$distance_to_snp)
	closest_gene <- unique(gene_info_copy[distance_to_snp==closest_gene_distance]$symbol)
	
	add <- data.table(snp=curr_snp, gene=closest_gene, p_value=curr_p, distance_to_snp=closest_gene_distance)
	res <- rbind(res,add)
	count <- count+1
}

res <- unique(res[!is.na(gene),])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write out
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite(res, sep="\t", "/Users/Liam/Desktop/nearest_gene_common_variants.tsv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Send to Minerva
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scp /Users/Liam/Desktop/nearest_gene_common_variants.tsv cottel02@chimera.hpc.mssm.edu:/sc/arion/projects/psychgen/methods/rx/data/
