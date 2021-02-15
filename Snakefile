import os
import pandas as pd
SampleTable = pd.read_table(config['sampletable'],index_col=0)
SAMPLES = list(SampleTable.index)
PAIRED_END= ('R2' in SampleTable.columns)
FRACTIONS= ['R1']
if PAIRED_END: FRACTIONS+= ['R2']

rule all:
    input:
        config["path"] + "output/seqs.fasta", 
        config["path"] + "output/seqtab.tsv", 
        config["path"] + "output/results.rds",
        config["path"] + "output/Nreads.tsv",
        config["path"] + "output/tax_gtdb.rds"

include: "rules/count.smk"        

rule QC:
    input:
        r1 = lambda wildcards: SampleTable.loc[wildcards.sample, 'R1'],
        r2 = lambda wildcards: SampleTable.loc[wildcards.sample, 'R2']
    output:
        r1 = config["path"] + "quality/{sample}_R1.fastq.gz",
        r2 = config["path"] + "quality/{sample}_R2.fastq.gz"
    params:
        adapters = os.path.abspath("../../Useful_Files/adapters.fa"),
    shell:
        "bbduk.sh in={input.r1} in2={input.r2} ref={params.adapters} out={output.r1} out2={output.r2} qtrim=r trimq=20  minlen=200 maq=20"

rule CutPrimers:
    input:
        r1 = config["path"] + "quality/{sample}_R1.fastq.gz",
        r2 = config["path"] + "quality/{sample}_R2.fastq.gz"
    output:
        r1 = config["path"] + "adapters/{sample}_R1.fastq.gz",
        r2 = config["path"] + "adapters/{sample}_R2.fastq.gz"
    params:
        primer_f = config["FORWARD"],
        primer_r = config["REVERSE"]
    shell:
        "cutadapt --discard-untrimmed -g ^{params.primer_f}  -G ^{params.primer_r} -o {output.r1} -p {output.r2} {input.r1} {input.r2}"   

        
rule dada2_filter:
    input:
        r1 =  expand(config["path"] + "adapters/{sample}_R1.fastq.gz",sample=SAMPLES),
        r2 =  expand(config["path"] + "adapters/{sample}_R2.fastq.gz",sample=SAMPLES)
    output:
        r1 = expand(config["path"] + "dada2_filtered/{sample}_R1.fastq.gz",sample=SAMPLES),
        r2 = expand(config["path"] + "dada2_filtered/{sample}_R2.fastq.gz",sample=SAMPLES),
        nreads= temp(config["path"] + "output/Nreads_filtered.txt")
    params:
        samples=SAMPLES
    threads:
        config['threads']
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/filter.txt"
    script:
        "scripts/dada2/filter.R"
        
rule learnErrorRates:
    input:
        r1= rules.dada2_filter.output.r1,
        r2= rules.dada2_filter.output.r2
    output:
        err_r1= config["path"] + "model/ErrorRates_r1.rds",
        err_r2 = config["path"] + "model/ErrorRates_r2.rds",
        plotErr1 = config["path"] + "figures/ErrorRates_r1.pdf",
        plotErr2 = config["path"] + "figures/ErrorRates_r2.pdf"
    threads:
        config['threads']
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/learnErrorRates.txt"
    script:
        "scripts/dada2/learnErrorRates.R"

rule dereplicate:
    input:
        r1 = rules.dada2_filter.output.r1,
        r2 = rules.dada2_filter.output.r2,
        err_r1 = rules.learnErrorRates.output.err_r1,
        err_r2 = rules.learnErrorRates.output.err_r2
    output:
        seqtab = temp(config["path"] + "output/seqtab_with_chimeras.rds"),
        nreads = temp(config["path"] + "output/Nreads_dereplicated.txt")
    params:
        samples = SAMPLES
    threads:
        config['threads']
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/dereplicate.txt"
    script:
        "scripts/dada2/dereplicate.R"

rule removeChimeras:
    input:
        seqtab = rules.dereplicate.output.seqtab
    output:
        seqtab = temp(config["path"] + "output/seqtab_nochimeras.rds"),
        nreads = temp(config["path"] + "output/Nreads_chimera_removed.txt")
    threads:
        config['threads']
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/removeChimeras.txt"
    script:
        "scripts/dada2/removeChimeras.R"

rule prepare_dbOTU:
    input:
        seqtab = rules.removeChimeras.output.seqtab
    output:
        fasta = temp(config["path"] + "output/seqs_dbOTU.fasta"),
        tsv = temp(config["path"] + "output/count_table.tsv")
    log:
        config["path"] + "logs/dada2/prepare_dbOTU.txt"
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dada2/dbOTU.R"
        
rule dbOTU:
    input:
        fasta = rules.prepare_dbOTU.output.fasta,
        tsv = rules.prepare_dbOTU.output.tsv
    output:
        tsv_out = temp(config["path"] + "output/dbOTU.tsv")
    log:
        log1 = config["path"] + "logs/dbOTU.log",
        log2 = config["path"] + "logs/dbOTU.debug"
    shell:
        "dbotu3.py --log {log.log1} --debug {log.log2} --output {output.tsv_out} {input.tsv} {input.fasta}"

rule import_dbOTU:
    input:
        dbOTU = rules.dbOTU.output.tsv_out,
        seqtab = rules.removeChimeras.output.seqtab
    output:
        rds = temp(config["path"] + "output/seqtab_dbOTU.rds")
    conda:
        "envs/dada2.yaml"
    log:
        log1 = config["path"] + "logs/import_dbOTU.log"
    script:
        "scripts/dada2/dbOTU_into.R"
        
        
rule filterLength:
    input:
        seqtab= rules.import_dbOTU.output.rds
    output:
        plot_seqlength = config["path"] + "figures/Lengths/Sequence_Length_distribution.pdf",
        plot_seqabundance = config["path"] + "figures/Lengths/Sequence_Length_distribution_abundance.pdf",
        rds= temp(config["path"] + "output/seqtab.rds"),
        tsv=  config["path"] + "output/seqtab.tsv",
    threads:
        1
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/filterLength.txt"
    script:
        "scripts/dada2/filterLength.R"

rule dada2_taxonomy:
    input:
        seqtab = rules.filterLength.output.rds
    output:
        tax = temp(config["path"] + "output/tax_silva.rds")
    params:
        silva = config['silva'],
        silva_species = config['silva_species']
    log:
        config["path"] + "logs/dada2/tax.txt"
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dada2/taxonomydada2.R"
        
rule ID_taxa:
    input:
        seqtab = rules.filterLength.output.rds,
    output:
        taxonomy= config["path"] + "output/tax_gtdb.rds"
    params:
        GTDB = config['GTDB']
    threads:
        config['threads']
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/IDtaxa_gtdb.txt"
    script:
        "scripts/dada2/IDtaxa.R"
        
        
rule dada2_end:
    input:
        seqtab = rules.filterLength.output.rds,
        tax = rules.dada2_taxonomy.output.tax
    output:
        fasta = config["path"] + "output/seqs.fasta",
        rds = config["path"] + "output/results.rds"
    conda:
        "envs/dada2.yaml"
    log:
        log1 = config["path"] + "logs/dada2_end.log"
    script:
        "scripts/dada2/dada2_end.R"
        
rule combine_read_counts:
    input:
        config["path"] + 'output/Nreads_raw.txt',
        config["path"] + 'output/Nreads_quality.txt',
        config["path"] + 'output/Nreads_adapters.txt',
        config["path"] + 'output/Nreads_filtered.txt',
        config["path"] + 'output/Nreads_dereplicated.txt',
        config["path"] + 'output/Nreads_chimera_removed.txt'
    output:
        config["path"] + 'output/Nreads.tsv',
    run:
        import pandas as pd
        import matplotlib
        import matplotlib.pylab as plt

        D = pd.read_table(input[0],sep=",",names=['00_raw'],index_col=0)
        D = D.join(pd.read_table(input[1],sep=",",names=['quality'],index_col=0))
        D = D.join(pd.read_table(input[2],sep=",",names=['adaptors'],index_col=0))
        D = D.join(pd.read_table(input[3],index_col=0))
        D = D.join(pd.read_table(input[4],index_col=0))
        D = D.join(pd.read_table(input[5],squeeze=True,index_col=0))
        D = D[['00_raw','quality','adaptors','filtered','denoised','merged','nonchim']]
        D.to_csv(output[0],sep='\t')
onsuccess:
    print(":) HAPPY")