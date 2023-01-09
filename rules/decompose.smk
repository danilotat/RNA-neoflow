rule decompose:
    conda: "../envs/vt.yaml"
    input:
        vcf_in = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"]+"/"+"{patient}.annotated.vcf"
    output:
        vcf_out = temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"]+"/"+"{patient}.decompose.vcf")
    log:
        config["OUTPUT_FOLDER"] + config["datadirs"]["logs"]["decompose"] + "/" + "{patient}.log"
    shell:
        "vt decompose -s {input.vcf_in} -o {output.vcf_out}"
