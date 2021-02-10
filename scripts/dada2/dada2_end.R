library(Biostrings)
library(tidyverse)

sink(snakemake@log[[1]])

seqtab <- readRDS(snakemake@input[['seqtab']]) 
tax <- readRDS(snakemake@input[['tax']])
seqtab %<>% t() %>% as.data.frame(stringsAsFactors=FALSE) %>% rownames_to_column(var="seqs")
tax %<>% as.data.frame(stringsAsFactors=FALSE) %>% rownames_to_column(var="seqs")%>% mutate(asv_id=paste0("asv",1:nrow(.)))
## Merge our data frame 
results <- left_join(tax,seqtab, by="seqs")
saveRDS(results, snakemake@output[['rds']])

fasta <- results$seqs; names(fasta) <- results$asv_id

writeXStringSet(DNAStringSet(fasta), snakemake@output[['fasta']],width=1000)