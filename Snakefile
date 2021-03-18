"""
Author: L. Velo Suarez
Affiliation: UMR-1078
Aim: A simple Snakemake workflow to process paired-end 16S amplicon from Illumina reads (MiSeq)
Run: snakemake --use-conda -j15   
Latest modification: 
  - todo
"""

configfile: "config.yaml"
report: "workflow.rst"

import os
import pandas as pd

SampleTable = pd.read_csv(config['sampletable'], sep='\t', index_col=False)
sample_id = list(SampleTable['sample_id'])

rule all:
    input:
        expand(config["path"] + "01_adapters/{sample}_R1.fastq.gz",sample=sample_id)
        #config["path"] + "QC/multiqc_report.html",
        #config["path"] + "output/results.fasta", 
        #config["path"] + "output/results.rds",
        #config["path"] + "output/seqtab_dbOTU.rds",
        #config["path"] + "output/tax_gtdb.rds",
        #config["path"] + "output/tax_silva.rds",
        #config["path"] + "output/Nreads.csv",
        
include: "rules/count.smk"
#include: "rules/qiime.smk"

def get_raw_fastq(sample_id):
    ''' return a dict with the path to the raw fastq files'''
    r1 = list(SampleTable[SampleTable["sample_id"] == sample_id]['R1'])[0].split(',')
    r2 = list(SampleTable[SampleTable["sample_id"] == sample_id]['R2'])[0].split(',')
    return {'r1': r1, 'r2': r2}

rule cutadapt:
    input: unpack(lambda wildcards: get_raw_fastq(wildcards.sample_id))
    output:
        fastq1 = config["path"] + "01_adapters/{sample_id}_R1.fastq.gz",
        fastq2 = config["path"] + "01_adapters/{sample_id}_R2.fastq.gz",
        qc = config["path"] + "QC/{sample_id}.qc.txt"
    params:
        adapters ="-g ^NCTACGGGNGGCWGCAG -G ^GACTACHVGGGGTATCTAATCC",
        extra = "--minimum-length 1 -q 20"
    log:
        config["path"] + "QC/{sample_id}_.log"
    wrapper:
        "0.72.0/bio/cutadapt/pe"

rule QC:
    input:
        r1 = config["path"] + "01_adapters/{sample}_R1.fastq.gz",
        r2 = config["path"] + "01_adapters/{sample}_R2.fastq.gz",
    output:
        r1 = config["path"] + "02_quality/{sample}_R1.fastq.gz",
        r2 = config["path"] + "02_quality/{sample}_R2.fastq.gz"
    params:
        adapters = os.path.abspath("../../Useful_Files/adapters.fa"),
    shell:
        "bbduk.sh in={input.r1} in2={input.r2} ref={params.adapters} out={output.r1} out2={output.r2} qtrim=rl trimq=20  minlen=200 maq=20"

rule fastqcR1:
    input: 
        config["path"] + "02_quality/{sample}_R1.fastq.gz"
    output: 
        html= config["path"] + "QC/{sample}_R1_fastqc.html",
        zip = config["path"] + "QC/{sample}_R1_fastqc.zip"
    wrapper: "0.72.0/bio/fastqc"

rule fastqcR2:
    input: 
        config["path"] + "02_quality/{sample}_R2.fastq.gz"
    output: 
        html= config["path"] + "QC/{sample}_R2_fastqc.html",
        zip = config["path"] + "QC/{sample}_R2_fastqc.zip"
    wrapper: "0.72.0/bio/fastqc"

rule multiqc:
    input: 
        expand([config["path"] + "QC/{sample}_R1_fastqc.html", 
                config["path"] + "QC/{sample}_R1_fastqc.zip",
                config["path"] + "QC/{sample}_R2_fastqc.html", 
                config["path"] + "QC/{sample}_R2_fastqc.zip",
                config["path"] + "QC/{sample}.qc.txt"],sample=sample_id)
    output: 
        config["path"] + "QC/multiqc_report.html"
    wrapper:  "0.72.0/bio/multiqc"
            
rule dada2_filter:
    input:
        r1 =  expand(config["path"] + "02_quality/{sample}_R1.fastq.gz",sample=sample_id),
        r2 =  expand(config["path"] + "02_quality/{sample}_R2.fastq.gz",sample=sample_id)
    output:
        r1 = expand(config["path"] + "03_dada2_filtered/{sample}_R1.fastq.gz",sample=sample_id),
        r2 = expand(config["path"] + "03_dada2_filtered/{sample}_R2.fastq.gz",sample=sample_id),
        nreads= temp(config["path"] + "output/Nreads_filtered.txt")
    params:
        samples=sample_id
    threads:
        config['threads']
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/filter.txt"
    script:
        "scripts/dada2/01_filter_dada.R"
        
rule learnErrorRates:
    input:
        r1 = expand(config["path"] + "03_dada2_filtered/{sample}_R1.fastq.gz",sample=sample_id),
        r2 = expand(config["path"] + "03_dada2_filtered/{sample}_R2.fastq.gz",sample=sample_id)
    output:
        err_r1= config["path"] + "04_model/ErrorRates_r1.rds",
        err_r2 = config["path"] + "04_model/ErrorRates_r2.rds",
        plotErr1 = report(config["path"] + "figures/ErrorRates_r1.png",category="QC reads"),
        plotErr2 = report(config["path"] + "figures/ErrorRates_r2.png",category="QC reads")
    threads:
        config['threads']
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/learnErrorRates.txt"
    script:
        "scripts/dada2/02_learnErrorRates.R"

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
        samples = sample_id
    threads:
        config['threads']
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/dereplicate.txt"
    script:
        "scripts/dada2/03_dereplicate.R"

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
        "scripts/dada2/04_removeChimeras.R"

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
        "scripts/dada2/05_dbOTU.R"
        
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
        rds = config["path"] + "output/seqtab_dbOTU.rds"
    conda:
        "envs/dada2.yaml"
    log:
        log1 = config["path"] + "logs/import_dbOTU.log"
    script:
        "scripts/dada2/06_dbOTU_into.R"
        
rule dada2_taxonomy:
    input:
        seqtab = rules.import_dbOTU.output.rds
    output:
        tax = config["path"] + "output/tax_silva.rds"
    params:
        silva = config['silva'],
        silva_species = config['silva_species']
    log:
        config["path"] + "logs/dada2/tax.txt"
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dada2/07_taxonomydada2.R"
        
rule ID_taxa:
    input:
        seqtab = rules.import_dbOTU.output.rds
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
        "scripts/dada2/08_IDtaxa.R"
                
rule filter_f:
    input:
        seqtab = rules.import_dbOTU.output.rds,
        tax = rules.dada2_taxonomy.output.tax
    output:
        plot_seqlength_nofilter = report(config["path"] + "figures/Sequence_Length_distribution.png", category="QC reads"),
        plot_seqabundance_nofilter = report(config["path"] + "figures/Sequence_Length_distribution_abundance.png", category ="QC reads"),
        plot_seqlength = report(config["path"] + "figures/Sequence_Length_distribution_filtered.png", category="QC reads"),
        plot_seqabundance = report(config["path"] + "figures/Sequence_Length_distribution_abundance_filtered.png", category ="QC reads"),
        fasta = config["path"] + "output/results.fasta",
        nreads= temp(config["path"] + "output/Nreads_length.txt"),
        rds = config["path"] + "output/results.rds",
    conda:
        "envs/dada2.yaml"
    log:
        config["path"] + "logs/dada2/filter.txt"
    script:
        "scripts/dada2/09_filter.R"
        
rule combine_read_counts:
    input:
        config["path"] + 'output/Nreads_raw.txt',
        config["path"] + 'output/Nreads_quality.txt',
        config["path"] + 'output/Nreads_adapters.txt',
        config["path"] + 'output/Nreads_filtered.txt',
        config["path"] + 'output/Nreads_dereplicated.txt',
        config["path"] + 'output/Nreads_chimera_removed.txt',
        config["path"] + "output/Nreads_length.txt"
    output:
        report(config["path"] + 'output/Nreads.csv', category="QC reads"),
        report(config["path"] + 'output/Nreads.html', category="QC reads")
    run:
        import pandas as pd
        import matplotlib
        import altair as alt

        D = pd.read_table(input[0],sep=",",names=['gatc'],index_col=0)
        D = D.join(pd.read_table(input[1],sep=",",names=['quality'],index_col=0))
        D = D.join(pd.read_table(input[2],sep=",",names=['adaptors'],index_col=0))
        D = D.join(pd.read_table(input[3],index_col=0))
        D = D.join(pd.read_table(input[4],index_col=0))
        D = D.join(pd.read_table(input[5],squeeze=True,index_col=0))
        D = D.join(pd.read_table(input[6],index_col=1))
        D = D[['gatc','quality','adaptors','filtered','denoised','merged','nonchim','tax_filter']]
        D.to_csv(output[0],sep=',')
        D['file_id'] = D.index.copy()
        D = pd.melt(D, id_vars=['file_id'], value_vars=['gatc','quality','adaptors','filtered','denoised','merged','nonchim','tax_filter'], var_name='steps', value_name='reads')
        input_dropdown = alt.binding_select(options=['gatc','quality','adaptors','filtered','denoised','merged','nonchim','tax_filter'])
        selection = alt.selection_single(fields=['steps'], bind=input_dropdown, name='QC reads')
        alt.Chart(D).mark_bar().encode(
            x=alt.X('file_id:N'),
            y='reads:Q',
            tooltip = [alt.Tooltip('reads:Q'),alt.Tooltip('file_id:N')]
        ).properties(width=1200,
                    height=300
        ).configure_axis(labelFontSize=15,
                         titleFontSize=20
        ).add_selection(selection).transform_filter(selection).save(output[1])


onsuccess:
    print(":) HAPPY")