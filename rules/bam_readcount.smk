# bad usage of wildcards here


rule bam_readcount:
    input:
        vcf_in = config["datadirs"]["VCF_out"]+"/"+"{patient}.decompose.vcf",
        bam_file = config["datadirs"]["BQSR_2"]+"/"+"{patient}_recal.pass2.bam"
    output:
        vcf_out = temp(config["datadirs"]["VCF_out"]+"/"+"{patient}_bam_readcount_indel.tsv"),
        vcf_out_2 = temp(config["datadirs"]["VCF_out"]+"/"+"{patient}_bam_readcount_snv.tsv")
    params:
        ref_fasta = config["resources"]["genome"],
        outdir = config["datadirs"]["VCF_out"],
        samp_name = "{patient}"
    conda:
        "../envs/bam_readcount.yaml"
    log:
        config["datadirs"]["logs"]["bam_readcount"] + "/" + "{patient}.log"
    shell:
        """
        python3 scripts/bam_readcount_helper.py \
        {input.vcf_in} {params.samp_name} {params.ref_fasta} {input.bam_file} {params.outdir}  
        """