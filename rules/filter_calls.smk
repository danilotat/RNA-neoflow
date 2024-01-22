import os

rule filtercalls:
    input:
        vcf = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"]+ '/' + "{patient}_workflow"+"/results/variants/variants.vcf.gz"
    output:
        vcf = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + "{patient}_DP_filt.vcf.gz",
        index = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + "{patient}_DP_filt.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yml"    
    threads: 
        config["params"]["samtools"]["threads"]
    shell:
        """
        bcftools view -e "GT='mis'" {input.vcf} |\
         bcftools view -i "FILTER='PASS' & (DP > 10) & (FORMAT/AD[0:0] > 3)" --threads {threads} |\
         bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule createTOML:
    input:
        config_main = configpath,
        toml_template = config["resources"]["vcfanno_toml"]
    output:
        toml_file = config["OUTPUT_FOLDER"] + "vcfanno.toml"
    conda:
        "../envs/cyvcf2.yml"
    shell:
        """
        python3 ../scripts/createTOML.py -y {input.config_main} - t {input.toml_template} -o {output.toml_file}
        """

rule vcfanno:
    input:
        toml_file = config["OUTPUT_FOLDER"] + "vcfanno.toml",
        vcf = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + '{patient}_DP_filt.vcf.gz'
    output:
        vcf = temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + '{patient}_annot.vcf.gz'),
        index = temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + '{patient}_annot.vcf.gz.tbi')
    params:
        vcfanno_binary = config["resources"]["vcfanno_binary"],
        extra = "--permissive-overlap"
    threads:
        config["params"]["vcfanno"]["threads"]
    conda:
        "../envs/samtool.yml"
    shell:
        """
        ./{params.vcfanno_binary} -p {threads} {params.extra} {input.toml_file} {input.vcf} |\
        bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule germProb:
    input:
        vcf = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + '{patient}_annot.vcf.gz',
        index = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + '{patient}_annot.vcf.gz.tbi'
    output:
        vcf = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + '{patient}_annot_germProb.vcf.gz'
    conda:
        "../envs/cyvcf2.yml"
    shell:
        """
        python3 ../scripts/germProb.py {input.vcf} {output.vcf}
        """

rule indexgermProb:
    input:
        vcf = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + '{patient}_annot_germProb.vcf.gz'
    output:
        index = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + '/' + '{patient}_annot_germProb.vcf.gz.tbi'
    conda:
        "../envs/samtool.yml"
    shell:
        """
        tabix -p vcf {input.vcf}
        """
    



