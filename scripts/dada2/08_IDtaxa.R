#!/usr/bin/Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
sink(snakemake@log[[1]])

seqtab= readRDS(snakemake@input[['seqtab']]) # seqtab
dna <- DNAStringSet(getSequences(seqtab)) # Create a DNAStringSet from the ASVs
load(snakemake@params[['GTDB']]) #trainingSet
ids <- IdTaxa(dna, trainingSet, strand="both", processors=snakemake@threads , verbose=FALSE) # use all processors

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab)
taxid %<>% as.data.frame() %>% rownames_to_column(var="seqs")
saveRDS(taxid,snakemake@output[['taxonomy']])