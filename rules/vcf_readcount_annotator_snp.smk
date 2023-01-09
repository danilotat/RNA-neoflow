rule vcf_readcount_annotator_snp:
    input:
        bam_snv = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"]+"/"+"{patient}_bam_readcount_snv.tsv",
        vcf_decomposed = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"]+"/"+"{patient}.decompose.vcf"
    output:
        vcf_out = temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"]+"/"+"{patient}_vcf_readcount_annotator_snv.vcf")
    params:
        dna_or_rna=config["params"]["vcf_readcount_annotator"]["dna_or_rna"],
        samp_name = "{patient}"
    conda:
        "../envs/vatools.yml"
    log:
        config["OUTPUT_FOLDER"] + config["datadirs"]["logs"]["vcf_readcount_annotator_snp"] + "/" + "{patient}.log"
    shell:
        "vcf-readcount-annotator {input.vcf_decomposed} {input.bam_snv} {params.dna_or_rna} \
        -s {params.samp_name} -t snv -o {output.vcf_out}"
