rule vcf_expression_annotator:
    input:
        vcf_in = config["datadirs"]["VCF_out"]+"/"+"{patient}_vcf_readcount_annotator_FINAL.vcf",
        expression_file = config["datadirs"]["expression"]+'/gene_expression.tsv'
    output:
        vcf_out = temp(config["datadirs"]["VCF_out"]+"/"+"{patient}_vcf_expression_annotator.vcf")
    conda:
        "../envs/vatools.yml"
    log:
        config["datadirs"]["logs"]["vcf_expression_annotator"] + "/" + "{patient}.log"
    shell:
        """
        vcf-expression-annotator -s {wildcards.patient} -i "genes" -e {wildcards.patient} -o {output.vcf_out} \
        {input.vcf_in} {input.expression_file} custom gene 
        """
