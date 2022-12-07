rule CombineSomaticAndGermline:
    input:
        vcf_in_somatic = config["datadirs"]["VCF_out"]+"/"+"{patient}.somatic.vcf",
        vcf_in_germline = config["datadirs"]["VCF_germ"]+"/"+"{patient}.filtered.vcf.gz",
        ref_fasta = ref_fasta
    output:
        vcf_out = config["datadirs"]["VCF_out"]+"/"+"{patient}_combined_somatic_plus_germline.vcf"
    params:
        samp_name = "{patient}",
        gatk3_path=config["resources"]["gatk3_jar_path"]
    shell:
        """
        java -jar {params.gatk3_path} -T CombineVariants -R {input.ref_fasta} \
        --variant {input.vcf_in_germline} --variant {input.vcf_in_somatic} \
        -o {output.vcf_out} --assumeIdenticalSamples
        """

rule SortCombinedVCF:
    input:
        vcf_in = config["datadirs"]["VCF_out"]+"/"+"{patient}_combined_somatic_plus_germline.vcf",
        reference_dict= ref_dict
    output:
        vcf_out = config["datadirs"]["VCF_out"]+"/"+"{patient}_sorted_combined_somatic_plus_germline.vcf"
    params:
        samp_name = "{patient}",
        picard_path=config["resources"]["picard_path"]
    shell:
        """
        java -jar {params.picard_path} SortVcf I={input.vcf_in} \
        O={output.vcf_out} SEQUENCE_DICTIONARY={input.reference_dict}
        """

rule PhaseVariants:
    input:
        vcf_in = config["datadirs"]["VCF_out"]+"/"+"{patient}_sorted_combined_somatic_plus_germline.vcf",
        ref_fasta = ref_fasta,
        bam_file = config["datadirs"]["BQSR_2"]+"/"+"{patient}_recal.pass2.bam"
    output:
        vcf_out = config["datadirs"]["VCF_out"]+"/"+"{patient}_phased.vcf"
    params:
        samp_name = "{patient}",
        gatk3_path=config["resources"]["gatk3_jar_path"]
    shell:
        """
        java -jar {params.gatk3_path} -T ReadBackedPhasing \
        -R {input.ref_fasta} -I {input.bam_file} --variant {input.vcf_in} \
        -L {input.vcf_in} -o {output.vcf_out}
        """

rule VepPhase:
    conda: "../envs/vep.yaml"
    input:
        vcf = config["datadirs"]["VCF_out"]+"/"+"{patient}_phased.vcf",
        cache = config["resources"]["vep_cache_dir"],
        plugins = config["resources"]["vep_plugin_dir"]
    output:
        vcfout = config["datadirs"]["VCF_out"]+"/"+"{patient}_annotated_phased.vcf.gz"
    params:
        assembly = config["params"]["annotate_variants"]["assembly"]
    shell:
        """ 
        vep --input_file {input.vcf} \
        --output_file {output.vcfout} \
        --format vcf --vcf --symbol --terms SO --tsl \
        --assembly {params.assembly} --offline --cache --dir_cache {input.cache} \
        --plugin Frameshift --plugin Wildtype --force_overwrite \
        --dir_plugins {input.plugins}
        """

rule tabix_vep:
    input:
        vcf_in = config["datadirs"]["VCF_out"]+"/"+"{patient}_annotated_phased.vcf.gz"
    output:
        vcf_out = config["datadirs"]["VCF_out"]+"/"+"{patient}_annotated_phased.vcf.tbi"
    conda:
        "../envs/tabix.yml"
    shell:
        """
        tabix -p vcf {input.vcf_in}
        """