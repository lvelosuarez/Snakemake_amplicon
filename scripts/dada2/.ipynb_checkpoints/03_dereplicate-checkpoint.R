#!/usr/bin/Rscript
suppressPackageStartupMessages(library(dada2))
sink(snakemake@log[[1]])
sam.names= snakemake@params[["samples"]]
filtFs = snakemake@input[['r1']]
filtRs = snakemake@input[['r2']]
load(snakemake@input[['err_r1']]) # errF
load(snakemake@input[['err_r2']]) #errR
names(filtFs) <- sam.names
names(filtRs) <- sam.names
derep_forward <- derepFastq(filtFs, verbose=TRUE)
derep_reverse <- derepFastq(filtRs, verbose=TRUE) 
names(derep_forward) <- sam.names
names(derep_reverse) <- sam.names

dada_forward <- dada(derep_forward, err=errF, pool="pseudo", multithread=TRUE, verbose=TRUE)
dada_reverse <- dada(derep_reverse, err=errR, pool="pseudo", multithread=TRUE, verbose=TRUE)
mergers <- mergePairs(dada_forward,derep_forward,dada_reverse,derep_reverse, verbose=TRUE)
## ---- seqtab ----
seqtab.all <- makeSequenceTable(mergers)
## ---- save seqtab ----
saveRDS(seqtab.all, snakemake@output[['seqtab']])
## get N reads
getNreads <- function(x) sum(getUniques(x))
track <- cbind( sapply(dada_forward, getNreads), sapply(mergers, getNreads))
colnames(track) <- c( "denoised", "merged")
write.table(track,snakemake@output[['nreads']],sep='\t')