#!/usr/bin/env Rscript
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
source("/uoa/scratch/shared/Soil_Microbiology_Group/libs/R/write_dataset.R")
source("/uoa/scratch/shared/Soil_Microbiology_Group/libs/R/phyloseq_taxonomy_imputation.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args[1] = "16S"
}

tax_path = "/uoa/scratch/shared/Soil_Microbiology_Group/Reference_databases/"

if(args[1]=="16S") {
  tax_file = paste0(tax_path, "16S/silva_nr99_v138.1_train_set.fa.gz")
} else if(args[1]=="AOA") {
  tax_file = paste0(tax_path, "AOA/amoA_AOA_tax_cluster.fasta")
} else if(args[1]=="AOB") {
  tax_file = paste0(tax_path, "AOB/amoA_AOB_tax_cluster.fasta")
} else if(args[1]=="COM") {
  tax_file = paste0(tax_path, "AOA/amoA_AOB_tax_cluster.fasta")
}

# Import seqtab created by dada2
seqtab.nochim <- readRDS(file = "./02_out/seqtab.Rds")
# Assign taxonomy to ASVs using reference database
taxa <- assignTaxonomy(seqtab.nochim, tax_file, multithread = TRUE)
# Import and match sample metadata
samdf <- read.csv("./00_ref/metadata.csv")
samdf <- samdf[order(samdf$sample),]
rownames(samdf) <- sort(rownames(seqtab.nochim))
# Merge sequence table, taxonomy and sample metadata into phyloseq object
ps <-
  phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows = FALSE),
    sample_data(samdf),
    tax_table(taxa)
  )
# Limit taxonomy table to 6 levels
tax_table(ps) <- subset(tax_table(ps), select = -c(8:12))
# Fix blank names and NA in tax table
ps <- phyloseq_taxonomy_imputation(ps)
# Add ASV sequences to phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
# Rename ASVs to include "ASV" as a prefix
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
# Write phyloseq object, ASV table, ASV fasta and ASV tax table to files
saveRDS(ps, file = paste0("./02_out/ps_", args[1], ".Rds"))
write.dataset(ps, filePATH = "./02_out/", filePREFIX = args[1])