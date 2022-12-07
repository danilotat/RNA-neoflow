rule salmon_quantification:
    input:
        unpack(get_fastq)
    output:
        quant = config["datadirs"]["salmon_quant"] + '/' + '{patient}' + '/quant.sf'
    params:
        index = config["resources"]["salmon_idx"],
        libtype = config["params"]["salmon"]["libtype"],
        zip_ext = config["params"]["salmon"]["zip_ext"],
        extra = config["params"]["salmon"]["extra"],
        outdir = config["datadirs"]["salmon_quant"] + '/' + '{patient}'
    threads: 
        config["params"]["threads"]["salmon"]
    conda:
        "../envs/salmon_new.yml"
    log:
        config["datadirs"]["logs"]["salmon_quant"] + '/{patient}.log'
    shell:
        """
        salmon quant -l {params.libtype} -i {params.index} -1 {input.r1} -2 {input.r2} -p {threads} -o {params.outdir}
        """

rule export_quantification:
    input:
        quant = expand(config["datadirs"]["salmon_quant"] + '/' + '{patient}' + '/quant.sf', patient=patients),
        index = config["resources"]["salmon_idx"],
        cdna_fasta = config["resources"]["transcriptome"],
        annotation = config["resources"]["gtf"]
    output:
        transcript = config["datadirs"]["expression"] + '/transcript_expression.tsv',
        gene = config["datadirs"]["expression"] + '/gene_expression.tsv'
    params:
        outdir = config["datadirs"]["expression"],
        patients = expand('{patient}', patient=patients)
    conda:
        "../envs/merge_salmon_quant.yml"
    script:
        "../scripts/merge_salmon_quantification.R"


