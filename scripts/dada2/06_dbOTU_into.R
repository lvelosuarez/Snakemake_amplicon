#!/usr/bin/Rscript
suppressPackageStartupMessages(library(tidyverse))
sink(snakemake@log[[1]])

seqtab= readRDS(snakemake@input[['seqtab']]) 

dbOTU <- seqtab %>% t() %>% as.data.frame(stringsAsFactors = FALSE) %>% rownames_to_column(var="seqs") %>% mutate(asv_id=paste0("asv", 1:nrow(.))) %>% dplyr::select(asv_id,seqs,everything())

table <- read_tsv(snakemake@input[['dbOTU']], col_names=TRUE)

dada2 <- data.frame(OTU_ID=dbOTU$asv_id,seq=dbOTU$seqs, stringsAsFactors= FALSE) %>% 
                          right_join(table, by="OTU_ID") %>% 
                          as.data.frame() %>% 
                          dplyr::select(-OTU_ID) %>% 
                          column_to_rownames(var = "seq") %>% 
                          t() %>%
                          as.matrix()
                
saveRDS(dada2, snakemake@output[['rds']])