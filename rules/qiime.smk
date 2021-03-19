# Snakemake rules imported in the main Snakefile to calculate SEPP trees for sequence results

rule import_fasta:
    input:
        fasta = config["path"] + "output/results.fasta"
    output:
        qiime_fasta = temp(config["path"] + "output/results.qza")
    conda:
        "envs/qiime2-2021.2-py36-linux-conda.yaml"
    shell:
        """ 
        qiime tools import --input-path {input.fasta} --output-path {output.qiime_fasta} --type 'FeatureData[Sequence]'
        """
rule SEPP:
    input:
        qiime_fasta = config["path"] + "output/results.qza"
    output:
        tree = temp(config["path"] + "output/insertion-tree.qza"),
        placements = temp(config["path"] + "output/insertion-placements.qza"),
    params:
        ref = os.path.abspath("../../Useful_Files/sepp-refs-gg-13-8.qza"),
    conda:
        "envs/qiime2-2021.2-py36-linux-conda.yaml"
    shell:
        """ 
        qiime fragment-insertion sepp --i-representative-sequences {input.qiime_fasta} --i-reference-database {params.ref} --o-tree {output.tree} --o-placements {output.placements}
        """
rule export_tree:
    input:
        tree = config["path"] + "output/insertion-tree.qza"),
    output:
        export = config["path"] + "output/tree.nwk"
    params:
        dir = config["path"]
    conda:
        "envs/qiime2-2021.2-py36-linux-conda.yaml"
    shell:
        """ 
        qiime tools export  --input-path {input.tree} --output-path {params.dir}
        """