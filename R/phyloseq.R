#!/usr/bin/env Rscript
library(dada2)
library(phyloseq, lib.loc = "/uoa/scratch/shared/Soil_Microbiology_Group/tools/r_libs/common")
library(Biostrings)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  args[1] = "16S"
}

tax_path = "/uoa/scratch/shared/Soil_Microbiology_Group/Reference_databases/dada2/"

if(args[1]=="AOA") {
  tax_file = paste0(tax_path, "amoA_AOA_tax_cluster.fasta")
} else if(args[1]=="AOB") {
  tax_file = paste0(tax_path, "amoA_AOB_tax_cluster.fasta")
} else if(args[1]=="16S") {
  tax_file = paste0(tax_path, "silva_nr99_v138.1_train_set.fa.gz")
} else {
  tax_file = paste0(tax_path, "amoA_AOA_tax_alves.fasta")
}

seqtab.nochim <- readRDS(file = "./02_out/seqtab.Rds")
taxa <- assignTaxonomy(seqtab.nochim, tax_file, multithread=TRUE)
samples.out <- rownames(seqtab.nochim)
samdf <- read.csv("./00_ref/metadata.csv")
rownames(samdf) <- samples.out
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
saveRDS(ps, file = "./02_out/ps.Rds")

OTU1 = t(as(otu_table(ps), "matrix"))
OTUdf = as.data.frame(OTU1)
write.table(OTUdf, file = "./02_out/ps.otu", sep = '\t')