# Moved from arcas-HLA to T1K, which is faster and runs smootly from conda through Snakemake.


rule genotype:
    input:
        unpack(get_fastq),
        idx=config["resources"]["t1k_file"],
    output:
        temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["HLA_typing"]
            + "/"
            + "{patient}_aligned_1.fa"
        ),
        temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["HLA_typing"]
            + "/"
            + "{patient}_aligned_2.fa"
        ),
        temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["HLA_typing"]
            + "/"
            + "{patient}_allele.tsv"
        ),
        temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["HLA_typing"]
            + "/"
            + "{patient}_allele.vcf"
        ),
        temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["HLA_typing"]
            + "/"
            + "{patient}_candidate_1.fq"
        ),
        temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["HLA_typing"]
            + "/"
            + "{patient}_candidate_2.fq"
        ),
        hla = config["OUTPUT_FOLDER"]
        + config["datadirs"]["HLA_typing"]
        + "/"
        + "{patient}_genotype.tsv",
    params:
        prefix="{patient}",
        outdir=lambda w, output: os.path.dirname(os.path.abspath(output.hla)),
    conda:
        "../envs/t1k.yml"
    threads: config["params"]["t1k"]["threads"]
    resources:
        time="4:00:00",
        ncpus=4,
        mem="32G",
    log:
        config["OUTPUT_FOLDER"] 
        + config["datadirs"]["logs"]["t1k"] 
        + "/" 
        + "{patient}.log",
    shell:
        """
        run-t1k -1 {input.r1} -2 {input.r2} --preset hla \
        -f {input.idx} -t {threads} -o {params.prefix} --od {params.outdir}
        """


rule extract_hla:
    input:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["HLA_typing"]
        + "/"
        + "{patient}_genotype.tsv",
    output:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["HLA_typing"]
        + "/"
        + "{patient}_allele_input_pvacseq.csv",
    conda:
        "../envs/cyvcf2.yml"
    log:
        config["OUTPUT_FOLDER"] 
        + config["datadirs"]["logs"]["t1k"]
        + "/"
        + "{patient}_hla.log",
    resources:
        time="0:20:00",
        ncpus=2,
        mem="8G",
    shell:
        "python3 scripts/HLA_typing.py {input} > {output}"
