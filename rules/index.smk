rule star_index:
    input: 
        fasta=config["resources"]["genome"],
        gtf=config["resources"]["gtf"]
    output:
        directory(config["datadirs"]["index_folder"]+"/genome_index")
    threads:
        config["params"]["thread"]
    conda:
        "../envs/star.yml"
    log:
        config["datadirs"]["logs"]["star_idx"] + '/{patient}.log'
    shell: "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} \
    --genomeFastaFiles {input.fasta} --sjdbOverhang 100 --sjdbGTFfile {input.gtf}"
