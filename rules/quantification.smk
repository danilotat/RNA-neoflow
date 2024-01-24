rule salmon_quantification:
    input:
        unpack(get_fastq)
    output:
        quant = config['OUTPUT_FOLDER'] + config["datadirs"]["salmon_quant"] + '/' + '{patient}' + '/quant.sf'
    params:
        index = config["resources"]["salmon_idx"],
        libtype = config["params"]["salmon"]["extra"]["libtype"],
        zip_ext = config["params"]["salmon"]["extra"]["zip_ext"],
        extra = config["params"]["salmon"]["extra"]["extra"],
        outdir = config['OUTPUT_FOLDER'] + config["datadirs"]["salmon_quant"] + '/' + '{patient}'
    threads: 
        config["params"]["salmon"]["threads"]
    resources:
        time="1:00:00",
        ncpus=4,
        mem="16G"
    conda:
        "../envs/salmon_new.yml"
    log:
        config['OUTPUT_FOLDER'] + config["datadirs"]["logs"]["salmon_quant"] + '/{patient}.log'
    shell:
        """
        salmon quant -l {params.libtype} -i {params.index} -1 {input.r1} -2 {input.r2} -p {threads} -o {params.outdir}
        """

rule export_quantification:
    input:
        quant = expand(config['OUTPUT_FOLDER'] + config["datadirs"]["salmon_quant"] + '/' + '{patient}' + '/quant.sf', patient=patients),
        index = config["resources"]["salmon_idx"],
        cdna_fasta = config["resources"]["transcriptome"],
        annotation = config["resources"]["gtf"]
    output:
        transcript = config['OUTPUT_FOLDER'] + config["datadirs"]["expression"] + '/transcript_expression.tsv',
        gene = config['OUTPUT_FOLDER'] + config["datadirs"]["expression"] + '/gene_expression.tsv'
    params:
        outdir = config['OUTPUT_FOLDER'] + config["datadirs"]["expression"],
        patients = expand('{patient}', patient=patients)
    conda:
        "../envs/merge_salmon_quant.yml"
    script:
        "../scripts/merge_salmon_quantification.R"


