#!/usr/bin/env Rscript
library(dada2)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args[1] = "16S"
  args[2] = "_R1"
  args[3] = "_R2"
}

trunc_len = 200
conc_status = TRUE

if(args[1]=="16S") {
  trunc_len_f = 220
  trunc_len_r = 220
} else if(args[1]=="AOA") {
  trunc_len_f = 200
  trunc_len_r = 200
  conc_status = TRUE
} else if(args[1]=="AOB") {
  trunc_len_f = 240
  trunc_len_r = 220
} else if(args[1]=="COM") {
  trunc_len_f = 200
  trunc_len_r = 200
}

inpath <- "./01_data/02_trimmed"
outpath <- "./01_data/03_filtered"

list.files(inpath)
fnFs <-
  sort(list.files(
    inpath,
    pattern = paste0(args[2], "_val_1.fq"),
    full.names = TRUE
  ))
fnRs <-
  sort(list.files(
    inpath,
    pattern = paste0(args[3], "_val_2.fq"),
    full.names = TRUE
  ))
if (length(fnFs) != length(fnRs))
  stop("Forward and reverse files do not match")
# Extract sample names using the forward read file name up to the forward read ID argument. eg SAMPLENAME_R1.fq.gz
sample.names <- sapply(strsplit(basename(fnFs), args[2]), `[`, 1)
# Place filtered files in filtered subdirectory
filtFs <-
  file.path(outpath, paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <-
  file.path(outpath, paste0(sample.names, "_2_filt.fastq.gz"))
# Name the filter objects by the sample names
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#Set truncLen and minLen according to your dataset. If too few reads are passing the filtering step, adjust maxEE up eg. c(2,5)
out <-
  filterAndTrim(
    fnFs,
    filtFs,
    fnRs,
    filtRs,
    maxN = 0,
    maxEE = c(2, 2),
    truncLen = c(trunc_len_f, trunc_len_r),
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
  )
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
# Plot error rates
errFplot <- plotErrors(errF, nominalQ = TRUE)
ggplot2::ggsave(errFplot, file = "02_out/filtF.error.rate.pdf")
errRplot <- plotErrors(errF, nominalQ = TRUE)
ggplot2::ggsave(errRplot, file = "02_out/filtR.error.rate.pdf")
# Infer sequence variants
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
# Merge read pairs
mergers <-
  mergePairs(dadaFs,
             filtFs,
             dadaRs,
             filtRs,
             justConcatenate = conc_status,
             verbose = TRUE)
seqtab <- makeSequenceTable(mergers)
# Remove chimeras. If too few reads are passing this step, check your adapter trimming
seqtab.nochim <-
  removeBimeraDenovo(seqtab,
                     method = "consensus",
                     multithread = TRUE,
                     verbose = TRUE)
saveRDS(seqtab.nochim, file = "./02_out/seqtab.Rds")
getN <- function(x)
  sum(getUniques(x))
track <-
  cbind(
    out,
    sapply(dadaFs, getN),
    sapply(dadaRs, getN),
    sapply(mergers, getN),
    rowSums(seqtab.nochim)
  )
# If processing a single sample, remove the sapply calls
#track <-
#  cbind(out,
#        getN(dadaFs),
#        getN(dadaRs),
#        getN(mergers),
#        rowSums(seqtab.nochim))
colnames(track) <-
  c("input",
    "filtered",
    "denoisedF",
    "denoisedR",
    "merged",
    "nonchim")
rownames(track) <- sample.names
write.csv(track, file = "./02_out/sample_tracking.csv")
