library(dada2)
sink(snakemake@log[[1]])

dada2= readRDS(snakemake@input[['seqtab']]) # seqtab.all

set.seed(100)
taxa <- assignTaxonomy(dada2,snakemake@params[['silva']] , multithread=TRUE, tryRC = TRUE,minBoot=50)
taxa.species <- addSpecies(taxa, snakemake@params[['silva_species']], verbose=TRUE)

saveRDS(taxa.species, snakemake@output[['tax']])