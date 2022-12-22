#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(dada2)
library(ggplot2)
inpath <- "./01_data/02_trimmed"
outpath <- "./01_data/03_filtered"
list.files(inpath)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_val_1.fastq and SAMPLENAME_R2_val_2.fq
fnFs <- sort(list.files(inpath, pattern = paste0(args[1], "_val_1.fq"), full.names = TRUE))
fnRs <- sort(list.files(inpath, pattern = paste0(args[2], "_val_2.fq"), full.names = TRUE))
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match")
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), args[1]), `[`, 1)
# Place filtered files in filtered subdirectory
filtFs <- file.path(outpath, paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(outpath, paste0(sample.names, "_2_filt.fastq.gz"))
# Name the filter objects by the sample names
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#Set truncLen and minLen according to your dataset
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,5), truncLen=c(240,160),truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# Plot error rates
errFplot <- plotErrors(errF, nominalQ=TRUE)
ggplot2::ggsave(errFplot, file = "02_out/filtF.error.rate.pdf")
errRplot <- plotErrors(errF, nominalQ=TRUE)
ggplot2::ggsave(errRplot, file = "02_out/filtR.error.rate.pdf")
# Infer sequence variants
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# Merge read pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, file = "./02_out/seqtab.Rds")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls
# track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file = "./02_out/sample_tracking.csv")
