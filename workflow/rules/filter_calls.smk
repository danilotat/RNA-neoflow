import os


rule vcfanno:
    input:
        toml_file=config["OUTPUT_FOLDER"] + "vcfanno.toml",
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_workflow"
        + "/results/variants/variants.vcf.gz",
    output:
        vcf=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["VCF_out"]
            + "/"
            + "{patient}_annot.vcf.gz"
        ),
        index=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["VCF_out"]
            + "/"
            + "{patient}_annot.vcf.gz.tbi"
        ),
    params:
        vcfanno_binary=config["resources"]["vcfanno_binary"],
        extra="--permissive-overlap",
        lua=config["resources"]["vcfanno_lua"],
    threads: config["params"]["vcfanno"]["threads"]
    resources:
        time="1:00:00",
        ncpus=4,
        mem="16G",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        {params.vcfanno_binary} -lua {params.lua} -p {threads} {params.extra} {input.toml_file} {input.vcf} |\
        bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule filtercalls:
    input:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot.vcf.gz",
    output:
        vcf=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["VCF_out"]
            + "/"
            + "{patient}_DP_filt.vcf.gz"
        ),
        index=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["VCF_out"]
            + "/"
            + "{patient}_DP_filt.vcf.gz.tbi"
        ),
    conda:
        "../envs/samtools.yml"
    threads: config["params"]["samtools"]["threads"]
    resources:
        time="0:20:00",
        ncpus=2,
        mem="8G",
    shell:
        """
        bcftools view -e "GT='mis'" {input.vcf} |\
         bcftools view -i "FILTER='PASS' & (DP > 5) & (FORMAT/AD[0:0] > 2)" --threads {threads} |\
         bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule createTOML:
    input:
        config_main=configpath,
        toml_template=config["resources"]["vcfanno_toml"],
    output:
        toml_file=config["OUTPUT_FOLDER"] + "vcfanno.toml",
    conda:
        "../envs/cyvcf2.yml"
    resources:
        time="0:20:00",
        ncpus=1,
        mem="1G",
    shell:
        """
        python3 scripts/createTOML.py -y {input.config_main} -t {input.toml_template} -o {output.toml_file}
        """


rule germProb:
    input:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_DP_filt.vcf.gz",
        index=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_DP_filt.vcf.gz.tbi",
    output:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot_germProb.vcf.gz",
    conda:
        "../envs/cyvcf2.yml"
    resources:
        time="0:20:00",
        ncpus=1,
        mem="4G",
    shell:
        """
        python3 scripts/germProb.py {input.vcf} {output.vcf}
        """


rule indexgermProb:
    input:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot_germProb.vcf.gz",
    output:
        index=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot_germProb.vcf.gz.tbi",
    conda:
        "../envs/samtools.yml"
    resources:
        time="0:20:00",
        ncpus=1,
        mem="4G",
    shell:
        """
        tabix -p vcf {input.vcf}
        """