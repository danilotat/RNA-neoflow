rule BQSR_1:
    input:
        bam=config["datadirs"]["bams"]+"/"+"{patient}_split.out.bam",
        GSNPs=config["resources"]["gsnps"],
        indel=config["resources"]["indel"],
        DbSNP=config["resources"]["dbsnps"],
        fasta=config["resources"]["genome"]
    output:
        recall=temp(config["datadirs"]["bams"] + "/" + "{patient}_recal.table")
    resources:
        mem_mb=config["params"]["RAM"]["BQSR"]
    threads: 
        config["params"]["threads"]["BQSR"]
    conda:
        "../envs/gatk.yml"
    log:
        config["datadirs"]["logs"]["base_recalibration"] + "/" + "{patient}.log"
    shell:
        """
        gatk BaseRecalibrator \
        -I {input.bam} \
        -R {input.fasta} \
        --known-sites  {input.GSNPs} \
        --known-sites  {input.indel}  \
        --known-sites  {input.DbSNP} \
        -O {output.recall}
        """

rule applyBQSR:
    input:
        bam=config["datadirs"]["bams"]+"/"+"{patient}_split.out.bam",
        fasta=config["resources"]["genome"],
        recall=config["datadirs"]["bams"] + "/" + "{patient}_recal.table"
    output:
        rbam=temp(config["datadirs"]["bams"]+'/'+"{patient}_recal.pass1.bam")
    threads: 
        config["params"]["threads"]["BQSR"]
    conda:
        "../envs/gatk.yml"
    resources:
        mem_mb=config["params"]["RAM"]["BQSR"]
    log:
        config["datadirs"]["logs"]["base_recalibration"] + "/" + "{patient}.log"
    shell:
        """
        gatk ApplyBQSR \
        -I {input.bam}  \
        -R {input.fasta} \
        --bqsr-recal-file {input.recall} \
        -O {output.rbam}
        """
rule BQSR_Pass2:         
    input:
        bam=config['datadirs']['bams'] + "/" + "{patient}_recal.pass1.bam",
        GSNPs=config["resources"]["gsnps"],
        indel=config["resources"]["indel"],
        DbSNP=config["resources"]["dbsnps"],
        fasta=config["resources"]["genome"]
    output:
        recall=temp(config['datadirs']['bams'] + "/" + "{patient}_recal_2.table")
    threads: 
        config["params"]["threads"]["BQSR"]
    conda:
        "../envs/gatk.yml"
    resources:
        mem_mb=config["params"]["RAM"]["BQSR"]
    log:
        config["datadirs"]["logs"]["base_recalibration"] + "/" + "{patient}.log"
    shell:
        """
        gatk BaseRecalibrator \
        -I {input.bam} \
        -R {input.fasta} \
        --known-sites  {input.GSNPs} \
        --known-sites  {input.indel}  \
        --known-sites  {input.DbSNP} \
        -O {output.recall}
        """ 

rule ApplyBQSR_2:
    # this is latest edited bam file, so we decide to keep it for further checks. 
    input:
        bam=config['datadirs']['bams'] + "/" + "{patient}_recal.pass1.bam",
        fasta=config["resources"]["genome"],
        recall=config['datadirs']['bams'] + "/" + "{patient}_recal_2.table"
    output:
        rbam=config['datadirs']['BQSR_2'] + "/" + "{patient}_recal.pass2.bam"
    threads: 
        config["params"]["threads"]["BQSR"]
    conda:
        "../envs/gatk.yml"
    resources:
        mem_mb=config["params"]["RAM"]["BQSR"]
    log:
        config["datadirs"]["logs"]["base_recalibration"] + "/" + "{patient}.log"
    shell:
        """
        gatk ApplyBQSR \
        -I {input.bam}  \
        -R {input.fasta} \
        --bqsr-recal-file {input.recall} \
        -O {output.rbam}
        """
