#!/usr/bin/Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ggplot2))
sink(snakemake@log[[1]])

set.seed(100)
errF <- learnErrors(snakemake@input[['r1']], nbases=snakemake@config[["learn_nbases"]], multithread=snakemake@threads,randomize = TRUE)
set.seed(100)
errR <- learnErrors(snakemake@input[['r2']], nbases=snakemake@config[["learn_nbases"]], multithread=snakemake@threads,randomize = TRUE)

save(errF,file=snakemake@output[['err_r1']])
save(errR,file=snakemake@output[['err_r2']])

## ---- plot-rates ----
plotErrors(errF,nominalQ=TRUE)
ggsave(snakemake@output[['plotErr1']])
plotErrors(errR,nominalQ=TRUE)
ggsave(snakemake@output[['plotErr2']])