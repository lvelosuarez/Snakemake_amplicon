# Snakemake rules imported in the main Snakefile to count read numbers in raw/quality/cutadapter R1 files.

rule count_raw_reads:
    input:
        r1 = SampleTable.R1.values
    output:
        nreads= temp(config["path"] + "output/Nreads_raw.txt")
    params:
        awk = """awk '!/file/{print $1,$4}'"""
    threads:
        15 #config['threads']
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | {params.awk} | sed 's/_R1.fastq.gz//' | sed 's/,//'| sed 's/ /,/' > {output.nreads}
        """
rule count_qc_reads:
    input:
        r1 = expand(config["path"] + "02_quality/{sample}_R1.fastq.gz",sample=sample_id)
    output:
        nreads= temp(config["path"] + "output/Nreads_quality.txt")
    params:
        awk = """awk '!/file/{print $1,$4}'"""
    threads:
        config['threads']
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | {params.awk} | sed 's/_R1.fastq.gz//' | sed 's/,//'| sed 's/ /,/' > {output.nreads}
        """
rule count_cutadapter_reads:
    input:
        r1 = expand(config["path"] + "01_adapters/{sample}_R1.fastq.gz",sample=sample_id)
    output:
        nreads= temp(config["path"] + "output/Nreads_adapters.txt")
    params:
        awk = """awk '!/file/{print $1,$4}'"""
    threads:
        config['threads']
    shell:
        """ 
        seqkit stats -j {threads} --basename {input.r1} | {params.awk}| sed 's/_R1.fastq.gz//' | sed 's/,//'| sed 's/ /,/' > {output.nreads}
        """
