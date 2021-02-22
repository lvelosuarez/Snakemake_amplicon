#!/usr/bin/Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))

sink(snakemake@log[[1]])

seqtab= readRDS(snakemake@input[['seqtab']]) # seqtab
tax <- readRDS(snakemake@input[['tax']])

seqtab %<>% t() %>% as.data.frame(stringsAsFactors=FALSE) %>% rownames_to_column(var="seqs")
tax %<>% as.data.frame(stringsAsFactors=FALSE) %>% rownames_to_column(var="seqs")%>% mutate(asv_id=paste0("asv",1:nrow(.)))

## Merge our data frame 
R <- left_join(tax,seqtab, by="seqs")
rm(seqtab,tax)
R_f <- R %>% filter(Kingdom == "Bacteria") %>% filter(!is.na(Phylum)) %>% filter(Family != "Mitochondria") 

# Length of sequences (how many times we can observe this length?)
## Get seq length number and plot (how many times we can observe this length?)
l_hist= as.data.frame(table(nchar(R$seqs)))
colnames(l_hist) <- c("length","count")
## plot l_hist
ggplot(l_hist,aes(x=length,y=count)) + 
  geom_col() + 
  ggtitle("Sequence Lengths by SEQ Count (no filter)") +
  theme_bw() +
    theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
## and save the graph
ggsave(snakemake@output[['plot_seqlength_nofilter']], width = 20, height = 8, units = "cm")
dev.off() 
## Get seq length abundance and plot (how many sequences have this length?)
a_hist <- tapply(rowSums(dplyr::select_if(R,is.numeric)), nchar(R$seqs),sum) %>% tibble("length"=names(.), "abundance"=.)
#plot a_hist
ggplot(a_hist,aes(x=length,y=abundance)) + 
  geom_col() + 
  ggtitle("Sequence Lengths by SEQ abundance (no filter)") +
  theme_bw() +
    theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
## and save the graph
ggsave(snakemake@output[['plot_seqabundance_nofilter']],, width = 20, height = 8, units = "cm")
dev.off() 
#################
#
# Filter by taxonomic assigments and replot
#
#################
## refaire graphs
lf_hist= as.data.frame(table(nchar(R_f$seqs)))
colnames(lf_hist) <- c("length","count")
## plot l_hist
ggplot(lf_hist,aes(x=length,y=count)) + 
  geom_col() + 
  ggtitle("Sequence Lengths by SEQ Count") +
  theme_bw() +
    theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
## and save the graph
ggsave(snakemake@output[['plot_seqlength']], width = 20, height = 8, units = "cm")
dev.off() 
## Get seq length abundance and plot (how many sequences have this length?)
b_hist <- tapply(rowSums(dplyr::select_if(R_f,is.numeric)), nchar(R_f$seqs),sum) %>% tibble("length"=names(.), "abundance"=.)
#plot b_hist
ggplot(b_hist,aes(x=length,y=abundance)) + 
  geom_col() + 
  ggtitle("Sequence Lengths by SEQ abundance") +
  theme_bw() +
    theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
## and save the graph
ggsave(snakemake@output[['plot_seqabundance']],, width = 20, height = 8, units = "cm")
dev.off() 
###########################
##
##
##
###########################
most_common_length <- as.numeric(b_hist[which.max(b_hist$abundance),'length'])
max_diff <- snakemake@config[['max_length_variation']]
right_length <- abs(nchar(R_f$seqs) - most_common_length) < max_diff
R_f <- R_f[right_length,]
##########################
saveRDS(R_f, snakemake@output[['rds']])

## Make a table with lost reads
table <- R %>% dplyr::select_if(is.numeric) %>% colSums() %>% tibble("samples"=names(.), "dbOTU"=.)
table <- R_f %>% dplyr::select_if(is.numeric) %>% colSums() %>% tibble("samples"=names(.), "tax_filter"=.) %>% left_join(table, by="samples") %>% dplyr::select(samples,dbOTU,tax_filter)
write.table(table,snakemake@output[['nreads']],sep='\t')
#### 
# Save seqs as fasta file for qiime
### 
fasta <- R_f$seqs; names(fasta) <- R_f$asv_id
writeXStringSet(DNAStringSet(fasta), snakemake@output[['fasta']],width=1000)