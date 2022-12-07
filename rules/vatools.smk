rule vatools:
    input:
        genin = config["datadirs"]["VCF_out"]+"/"+"{patient}.vcf"
    params: 
        samp_name = "{patient}"
    output:
        genout = temp(config["datadirs"]["VCF_out"]+"/"+"{patient}.genotype.vcf")
    conda:
        "../envs/vatools.yml"
    log:
        config["datadirs"]["logs"]["vatools"] + "/" + "{patient}.log"
    shell:
        "vcf-genotype-annotator {input.genin} {params.samp_name} 0/1 -o {output.genout}"
