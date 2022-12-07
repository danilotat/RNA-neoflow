rule HaplotypeCaller:
    input:
        bam=config["datadirs"]["BQSR_2"]+'/'+'{patient}_recal.pass2.bam',
        intervals=config["datadirs"]["utils"]+"/interval-files/"+"{interval}-scattered.interval_list",
        ref_fasta=ref_fasta,
        dbSNP=config["resources"]["dbsnps"]
    output:
        vcf=temp(config["datadirs"]["VCF_germ"]+'/{patient}.{interval}.unfiltered.vcf.gz'),
        vcf_tbi=temp(config["datadirs"]["VCF_germ"]+'/{patient}.{interval}.unfiltered.vcf.gz.tbi')
    params:
        threshold=config["params"]["HaplotypeCaller"]["threshold"]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk HaplotypeCaller -R {input.ref_fasta} \
        -L {input.intervals} -I {input.bam} \
        -O {output.vcf} \
        --dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling {params.threshold} \
        --dbsnp {input.dbSNP}
        """

rule merge_germ:
    input:
        get_mergevcfs_germ_input
    output:
        vcf=config["datadirs"]["VCF_germ"]+'/{patient}.unfiltered.vcf.gz'
    params:
        i=lambda wildcards, input: ['-I ' + vcf for vcf in input]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk MergeVcfs {params.i} -O {output.vcf}
        """

rule tabix_merged:
    input:
        vcf=config["datadirs"]["VCF_germ"]+'/{patient}.unfiltered.vcf.gz'
    output:
        vcf=config["datadirs"]["VCF_germ"]+'/{patient}.unfiltered.vcf.gz.tbi'
    conda:
        "../envs/tabix.yml"
    shell:
        """
        tabix -p vcf {input.vcf_in}
        """

rule VariantFiltration:
    input:
        ref_fasta=ref_fasta,
        vcf=config["datadirs"]["VCF_germ"]+'/{patient}.unfiltered.vcf.gz'
    output:
        filtered_vcf=config["datadirs"]["VCF_germ"]+'/{patient}.filtered.vcf.gz',
        filtered_vcf_idx=config["datadirs"]["VCF_germ"]+'/{patient}.filtered.vcf.gz.tbi'
    params:
        window=config["params"]["VariantFiltration"]["window"],
        cluster=config["params"]["VariantFiltration"]["cluster"],
        FS=config["params"]["VariantFiltration"]["FS"],
        QD=config["params"]["VariantFiltration"]["QD"]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk VariantFiltration --R {input.ref_fasta} \
        --V {input.vcf} --window {params.window} \
        --cluster {params.cluster} \
        --filter-name "FS" --filter {params.FS} \
        --filter-name "QD" --filter {params.QD} \
        -O {output.filtered_vcf} 
        """

    


    