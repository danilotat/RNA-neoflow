# Moved from arcas-HLA to T1K, which is faster and runs smootly from conda through Snakemake.

rule genotype:
    input:
        unpack(get_fastq),
        idx=config["datadirs"]["utils"]+"/"+"hlaidx_rna_seq.fa"
    output:
        temp(config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_aligned_1.fa'),
        temp(config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_aligned_2.fa'),
        temp(config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_allele.tsv'),
        temp(config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_allele.vcf'),
        temp(config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_candidate_1.fq'),
        temp(config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_candidate_2.fq'),
        config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_genotype.tsv'
    params:
        prefix='{patient}',
        outdir=config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']
    conda:
        "../envs/t1k.yml"
    threads: config['params']['threads']['t1k']
    log: config["OUTPUT_FOLDER"] + config["datadirs"]['logs']['t1k']
    shell:
        """
        run-t1k -1 {input.r1} -2 {input.r2} --preset hla \
        -f {input.idx} -t {threads} -o {params.prefix} --od {params.outdir}
        """



rule extract_hla:
    input:
        config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_genotype.tsv'
    output:
        config["OUTPUT_FOLDER"] + config["datadirs"]['HLA_typing']+'/'+'{patient}_allele_input_pvacseq.csv'
    shell:
        "python3 ../scripts/HLA_typing {input} > {output}"

