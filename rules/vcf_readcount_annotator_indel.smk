rule vcf_readcount_annotator_indel:
    input:
        bam_indel = config["datadirs"]["VCF_out"]+"/"+"{patient}_bam_readcount_indel.tsv",
        vcf_in = config["datadirs"]["VCF_out"]+"/"+"{patient}_vcf_readcount_annotator_snv.vcf"
    output:
        vcf_out = temp(config["datadirs"]["VCF_out"]+"/"+"{patient}_vcf_readcount_annotator_FINAL.vcf")
    params:
        dna_or_rna=config["params"]["vcf_readcount_annotator"]["dna_or_rna"],
        samp_name = "{patient}"
    conda:
        "../envs/vatools.yml"
    log:
        config["datadirs"]["logs"]["vcf_readcount_annotator_indel"] + "/" + "{patient}.log"
    shell:
        "vcf-readcount-annotator {input.vcf_in} {input.bam_indel} {params.dna_or_rna} \
        -s {params.samp_name} -t indel -o {output.vcf_out}"
