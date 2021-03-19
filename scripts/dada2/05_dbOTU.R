#!/usr/bin/Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))

sink(snakemake@log[[1]])

seqtab= readRDS(snakemake@input[['seqtab']]) 
dbOTU <- seqtab %>% t() %>% as.data.frame(stringsAsFactors = FALSE) %>% rownames_to_column(var="seqs") %>% mutate(asv_id=paste0("asv", 1:nrow(.))) %>% dplyr::select(asv_id,seqs,everything())

fasta <- dbOTU$seqs; names(fasta) <- dbOTU$asv_id

writeXStringSet(DNAStringSet(fasta), snakemake@output[['fasta']],width=1000)
write_tsv(dplyr::select(dbOTU,-seqs), snakemake@output[['tsv']])