rule annotate_variants:
    conda:
        "../envs/vep.yml"
    input:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_overlap.vcf.gz",
        #cache=config["resources"]["vep_cache_dir"],
        #plugins=config["resources"]["vep_plugin_dir"],
    output:
        vcfout=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}.annotated.vcf",
    params:
        assembly=config["params"]["annotate_variants"]["assembly"],
        filtering=config["params"]["annotate_variants"]["filtering"],
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["annotate_variants"]
        + "/"
        + "{patient}.log",
    shell:
        """ 
        vep --input_file {input.vcf} \
        --output_file {output.vcfout} \
        --format vcf --vcf --symbol --terms SO --tsl \
        -{params.filtering} \
        --plugin Frameshift --plugin Wildtype --force_overwrite \
        --assembly {params.assembly}
        """
        # --assembly {params.assembly} --offline --cache --dir_cache {input.cache} \
        # --dir_plugins {input.plugins}
        
