rule align:
    input:
        unpack(get_fastq),
        index = config["resources"]["star_index"]
    output:
        temp(config['OUTPUT_FOLDER'] + config["datadirs"]["mapped_reads"] + "/" + "{patient}_Aligned.sortedByCoord.out.bam")
        # temp(dir(config["datadirs"]["mapped_reads"] + "/" + "{patient}__STARgenome")),
        # temp(dir(config["datadirs"]["mapped_reads"] + "/" + "{patient}__STARpass1"))
    conda:
        "../envs/star.yml"
    params:
        index = lambda wc, input: input.index,
        prefix = config['OUTPUT_FOLDER'] + config["datadirs"]["mapped_reads"] + '/{patient}_',
        extra = "--sjdbGTFfile {} {}".format(config["resources"]["gtf"], config["params"]["STAR"]["extra"])
    threads:
        config["params"]["STAR"]["threads"]
    resources:
        mem = "40G",
        time = "10:00:00",
        ncpus = 8
    log:
        config['OUTPUT_FOLDER'] + config["datadirs"]["logs"]["align"] + "/" + "{patient}.log"
    shell:
        '''
        STAR --readFilesIn {input.r1} {input.r2} \
        --genomeDir {input.index} --runThreadN {threads} \
        --outFileNamePrefix {params.prefix} {params.extra}
        '''
