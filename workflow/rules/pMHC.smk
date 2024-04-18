import os

rule pMHCpeptides:
    input:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot_germProb.vcf.gz",
        vcf_idx=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot_germProb.vcf.gz.tbi",
        hla=config["OUTPUT_FOLDER"]
        + config["datadirs"]["HLA_typing"]
        + "/"
        + "{patient}_genotype.tsv",
    params:
        patname="{patient}",
        outfolder=lambda w, output: os.path.dirname(os.path.abspath(output.out)),
        threads=config["params"]["pMHC"]["threads"],
    output:
        out=config["OUTPUT_FOLDER"]
        + config["datadirs"]["peptides"]
        + "/"
        + "{patient}.epitopes.csv",
    container:
        "docker://danilotat/netmhcpan-minimal"
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["pMHC"]
        + "/"
        + "{patient}.log",
    resources:
        time="2:00:00",
        ncpus=4,
        mem="8G",
    shell:
        """
        python3 netmhcpan_launcher.py --vcf {input.vcf} -a {input.hla} -p {params.patname} -o {params.outfolder} -t {params.threads}
        """
