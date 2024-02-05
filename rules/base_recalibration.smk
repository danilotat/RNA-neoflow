rule BQSR_1:
    input:
        bam=config['OUTPUT_FOLDER'] + config["datadirs"]["bams"]+"/"+"{patient}_split.out.bam",
        GSNPs=config["resources"]["gsnps"],
        indel=config["resources"]["indel"],
        DbSNP=config["resources"]["dbsnp146"],
        fasta=config["resources"]["genome"]
    output:
        recall=config['OUTPUT_FOLDER'] + config["datadirs"]["bams"] + "/" + "{patient}_recal.table"
    resources:
        time="6:00:00",
        ncpus=4,
        mem="32G"
    threads: 
        config["params"]["BQSR"]["threads"]
    conda:
        "../envs/gatk.yml"
    log:
        config['OUTPUT_FOLDER'] + config["datadirs"]["logs"]["base_recalibration"] + "/" + "{patient}.log"
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
        bam=config['OUTPUT_FOLDER'] + config["datadirs"]["bams"]+"/"+"{patient}_split.out.bam",
        fasta=config["resources"]["genome"],
        recall=config['OUTPUT_FOLDER'] + config["datadirs"]["bams"] + "/" + "{patient}_recal.table"
    output:
        rbam=config['OUTPUT_FOLDER'] + config["datadirs"]["BQSR"]+'/'+"{patient}_recal.bam"
    threads: 
        config["params"]["BQSR"]["threads"]
    conda:
        "../envs/gatk.yml"
    resources:
        time="6:00:00",
        ncpus=4,
        mem="32G"
    log:
        config['OUTPUT_FOLDER'] + config["datadirs"]["logs"]["base_recalibration"] + "/" + "{patient}.log"
    shell:
        """
        gatk ApplyBQSR \
        -I {input.bam}  \
        -R {input.fasta} \
        --bqsr-recal-file {input.recall} \
        -O {output.rbam}
        """

### --- DO WE NEED 2 STEP OF BQSR ?? --- ### 


# rule BQSR_Pass2:         
#     input:
#         bam=config["OUTPUT_FOLDER"] + config["datadirs"]['bams'] + "/" + "{patient}_recal.pass1.bam",
#         GSNPs=config["resources"]["gsnps"],
#         indel=config["resources"]["indel"],
#         DbSNP=config["resources"]["dbsnps"],
#         fasta=config["resources"]["genome"]
#     output:
#         recall=temp(config["OUTPUT_FOLDER"] + config["datadirs"]['bams'] + "/" + "{patient}_recal_2.table")
#     threads: 
#         config["params"]["BQSR"]["threads"]
#     conda:
#         "../envs/gatk.yml"
#     resources:
#         mem_mb=config["params"]["BQSR"]["RAM"]
#     log:
#         config['OUTPUT_FOLDER'] + config["datadirs"]["logs"]["base_recalibration"] + "/" + "{patient}.log"
#     shell:
#         """
#         gatk BaseRecalibrator \
#         -I {input.bam} \
#         -R {input.fasta} \
#         --known-sites  {input.GSNPs} \
#         --known-sites  {input.indel}  \
#         --known-sites  {input.DbSNP} \
#         -O {output.recall}
#         """ 

# rule ApplyBQSR_2:
#     # this is latest edited bam file, so we decide to keep it for further checks. 
#     input:
#         bam=config["OUTPUT_FOLDER"] + config["datadirs"]['bams'] + "/" + "{patient}_recal.pass1.bam",
#         fasta=config["resources"]["genome"],
#         recall=config["OUTPUT_FOLDER"] + config["datadirs"]['bams'] + "/" + "{patient}_recal_2.table"
#     output:
#         rbam=config["OUTPUT_FOLDER"] + config["datadirs"]['BQSR_2'] + "/" + "{patient}_recal.pass2.bam"
#     threads: 
#         config["params"]["BQSR"]["threads"]
#     conda:
#         "../envs/gatk.yml"
#     resources:
#         mem_mb=config["params"]["BQSR"]["RAM"]
#     log:
#         config['OUTPUT_FOLDER'] + config["datadirs"]["logs"]["base_recalibration"] + "/" + "{patient}.log"
#     shell:
#         """
#         gatk ApplyBQSR \
#         -I {input.bam}  \
#         -R {input.fasta} \
#         --bqsr-recal-file {input.recall} \
#         -O {output.rbam}
#         """
