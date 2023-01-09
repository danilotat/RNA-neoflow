rule ref_transcript_mismatch_reporter:
    input:
        vcf_in = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_vcf_expression_annotator.vcf" 
    params:
        filter_vcf = config["params"]["ref_transcript_mismatch_reporter"]
    output:
        vcf_out = temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_ref_transcript_mismatch_reporter.vcf")
    conda:
        "../envs/vatools.yml"
    log:
        config["OUTPUT_FOLDER"] + config["datadirs"]["logs"]["ref_transcript_mismatch_reporter"]+'/'+"{patient}.log"
    shell:
        "ref-transcript-mismatch-reporter {input.vcf_in} --filter {params.filter_vcf} -o {output.vcf_out}"

rule bgzip_mismatch_reporter:
    input:
        vcf_in = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_ref_transcript_mismatch_reporter.vcf"
    output:
        vcf_out = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_ref_transcript_mismatch_reporter.vcf.gz"
    conda:
        "../envs/samtool.yml"
    shell:
        """
        bgzip -c {input.vcf_in} > {output.vcf_out}
        """
rule tabix_mismatch_reporter:
    input:
        vcf_in = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_ref_transcript_mismatch_reporter.vcf.gz"
    output:
        vcf_out = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_ref_transcript_mismatch_reporter.vcf.gz.tbi"
    conda:
        "../envs/tabix.yml"
    shell:
        """
        tabix -p vcf {input.vcf_in}
        """




