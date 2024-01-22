rule Strelka_prep:
    input:
        bam=config['OUTPUT_FOLDER'] + config["datadirs"]["BQSR"]+'/'+"{patient}_recal.bam"
    output:
        config["OUTPUT_FOLDER"] + config["datadirs"]['VCF_out']+'/'+'{patient}_workflow'+'/'+'runWorkflow.py'
    params:
        ref_fasta=config['resources']['genome'],
        regions=config['resources']['intervals_coding'],
        runDir=config["OUTPUT_FOLDER"] + config["datadirs"]['VCF_out']+'/'+'{patient}_workflow'
    conda:
        "../envs/strelka2.yml"
    shell:
        """
        configureStrelkaGermlineWorkflow.py \
        --bam {input.bam} \
        --rna \
        --referenceFasta {params.ref_fasta} \
        --callRegions {params.regions} \
        --runDir {params.runDir} \
        --reportEVSFeatures
        """
rule Strelka2:
    input:
        script=config["OUTPUT_FOLDER"] + config["datadirs"]['VCF_out']+'/'+'{patient}_workflow'+'/'+'runWorkflow.py'
    output:
        config["OUTPUT_FOLDER"] + config["datadirs"]['VCF_out']+'/'+'{patient}_workflow'+'/results/variants/'+'variants.vcf.gz'
    params:
        threads=config['params']['strelka2']['threads']
    conda:
        "../envs/strelka2.yml"
    shell:
        """
        {input.script} -m local -j {params.threads}
        """
