rule sort_bam:
    input:
        bam_in = config["datadirs"]["mapped_reads"] + "/" + "{patient}_Aligned.out.bam"
    output:
        bam_out = temp(config["datadirs"]["mapped_reads"] + "/" + "{patient}_Aligned.sortedByCoord.out.bam")
    conda:
        "../envs/samtool.yml"
    log:
        config["datadirs"]["logs"]["sort_bam"] + "/" + "{patient}.log"
    shell:
        '''
        samtools sort {input.bam_in} -o {output.bam_out}
        '''