rule star_index:
    input:
        fasta=config["resources"]["genome"],
        gtf=config["resources"]["gtf"],
    output:
        directory(config["datadirs"]["index_folder"] + "/genome_index"),
    threads: config["params"]["thread"]
    conda:
        "../envs/star.yml"
    log:
        config["datadirs"]["logs"]["star_idx"] + "/{patient}.log",
    resources:
        mem="60G",
        ncpus=8,
        time="6:00:00",
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} \
    --genomeFastaFiles {input.fasta} --sjdbOverhang 100 --sjdbGTFfile {input.gtf}"


rule salmon_gentrome:
    input:
        genome=config["resources"]["genome"],
        cdna=config["resources"]["transcriptome"],
    output:
        temp(config["resources"] + "/" + "gentrome.fa.gz"),
    shell:
        "cat {input.genome} {input.cdna} | gzip > {output}"


rule salmon_idx:
    input:
        gentrome=config["resources"] + "/" + "gentrome.fa.gz",
    output:
        config["resources"]["salmon_idx"] + "/" + "ctable.bin",
    threads: config["params"]["threads"]["salmon"]
    params:
        outdir=config["resources"]["salmon_idx"],
        extra=config["params"]["salmon"]["index"],
    resources:
        mem="40G",
        ncpus=8,
        time="4:00:00",
    conda:
        "../envs/salmon_new.yml"
    log:
        config["datadirs"]["logs"]["salmon_quant"] + "/" + "index.log",
    shell:
        """
        salmon index -t {input.gentrome} -i {params.outdir} \
        --decoys {input.decoys} {params.extra}
        """
