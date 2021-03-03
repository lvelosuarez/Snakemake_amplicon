#!/usr/bin/Rscript
suppressPackageStartupMessages(library(dada2))
sink(snakemake@log[[1]])

seqtab.all= readRDS(snakemake@input[['seqtab']]) # seqtab.all


# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab.all, method=snakemake@config[['chimera_method']], multithread=snakemake@threads,verbose=TRUE)
seqtable_nochim_collapse <- collapseNoMismatch(seqtab, verbose=TRUE)
saveRDS(seqtable_nochim_collapse, snakemake@output[['seqtab']])

track <- rowSums(seqtable_nochim_collapse)
names(track) <- row.names(seqtable_nochim_collapse)

write.table(track,col.names = c("nonchim"),
            snakemake@output[['nreads']],sep='\t')