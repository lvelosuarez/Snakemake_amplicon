# Snakemake rules imported in the main Snakefile to filter Ns in R2 files of gatc lib1 (february 2021)



if config['option'] == 'LIB2':
    rule two_a:
        input:
            'foo.txt',
        output:
            'bar.txt',
elif config['option'] == 'LIB1':
    rule two_b:
        input:
            'foo.txt',
        output:
            'bar.txt',
else:
    sys.exit() ## handle this case