library(dada2)

sink(snakemake@log[[1]])

track.filt <- filterAndTrim(snakemake@input[['r1']],snakemake@output[['r1']], 
                            snakemake@input[['r2']],snakemake@output[['r2']], 
                            maxEE=snakemake@config[["maxEE"]], 
                            compress=TRUE,
                            verbose=TRUE,
                            multithread=snakemake@threads)
row.names(track.filt) <- snakemake@params[["samples"]]
colnames(track.filt) = c('raw','filtered')
write.table(track.filt,snakemake@output[['nreads']],  sep='\t')