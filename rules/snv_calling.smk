rule Mutect2:
    input:
        reference=config["resources"]["mutect_reference"],
        genome=config["resources"]["genome"],
        tumor=config["OUTPUT_FOLDER"] + config["datadirs"]["BQSR_2"]+'/'+'{patient}_recal.pass2.bam',
        intervals=config["OUTPUT_FOLDER"] + config["datadirs"]["utils"]+"/interval-files/"+"{interval}-scattered.interval_list",
        PoN=config["resources"]["PoN"]
    params:
        extra=config["params"]["mutect2"]
    output:
        vcf=temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.{interval}.unfiltered.vcf.gz"),
        f1r2=temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.{interval}.f1r2.tar.gz"),
        stats=temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.{interval}.unfiltered.vcf.gz.stats")
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk Mutect2 -R {input.genome} \
        -I {input.tumor} \
        --germline-resource {input.reference} \
        {params.extra} \
        -L {input.intervals} \
        -pon {input.PoN} \
        --f1r2-tar-gz {output.f1r2} -O {output.vcf}
        """

rule pileup_summaries:
    input:
        bam=config["OUTPUT_FOLDER"] + config["datadirs"]["BQSR_2"]+'/'+'{patient}_recal.pass2.bam',
        germ_res=config["resources"]["contamination_resource"]
    output:
        pileup_resume=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}_pileupsummaries.table"
    conda:
        "../envs/gatk.yml"
    log:
        config["OUTPUT_FOLDER"] + config["datadirs"]["logs"]["snv_calling"]+'/'+"{patient}.log"
    shell:
        """
        gatk GetPileupSummaries -I {input.bam} -V {input.germ_res} -L {input.germ_res} \
        -O {output.pileup_resume}
        """

rule calculate_contamination:
    input:
        pileup_resume=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}_pileupsummaries.table"
    output:
        table=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}_contamination.table"
    conda:
        "../envs/gatk.yml"
    log:
        config["OUTPUT_FOLDER"] + config["datadirs"]["logs"]["snv_calling"]+'/'+"{patient}.log"
    shell:
        "gatk CalculateContamination -I {input.pileup_resume} -O {output.table}"

rule merge_vcfs:
    input:
        get_mergevcfs_input
    output:
        vcf=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.unfiltered.vcf.gz"
    params:
        i=lambda wildcards, input: ['-I ' + vcf for vcf in input]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk MergeVcfs {params.i} -O {output.vcf}
        """

rule merge_stats:
    input:
        get_mergestats_input
    output:
        stats=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}_unfiltered.vcf.stats"
    params:
        i=lambda wildcards, input: ['-stats ' + s for s in input]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk MergeMutectStats {params.i} -O {output.stats} 
        """

rule filtering_calls:
    input:
        vcf_unfiltered=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.unfiltered.vcf.gz",
        ref=config["resources"]["genome"],
        contamination=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}_contamination.table",
        stats=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}_unfiltered.vcf.stats"
        # f1r2_model=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}_read_orientation_model.tar.gz"
    output:
        filtered_vcf=temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.filtered.vcf"),
        idx=temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.filtered.vcf.idx"),
        stats=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.filtered.vcf.filteringStats.tsv"
    params:
        extra=config["params"]["FilterMutectCalls"]
    conda:
        "../envs/gatk.yml"
    log:
        config["OUTPUT_FOLDER"] + config["datadirs"]["logs"]["snv_calling"]+'/'+"{patient}.log"
    shell:
        """
        gatk FilterMutectCalls -V {input.vcf_unfiltered} -R {input.ref} \
        --contamination-table {input.contamination} --stats {input.stats} \
        -O {output.filtered_vcf} \
        {params.extra}
        """

rule tabix_filtered_calls:
    input:
        vcf=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.filtered.vcf"
    output:
        vcf_gz=temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.filtered.vcf.gz"),
        vcf_idx=temp(config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]+"/{patient}.filtered.vcf.gz.tbi")
    conda:
        '../envs/tabix.yml'
    shell:
        '''
        bgzip {input.vcf}
        tabix -p vcf {output.vcf_gz}
        '''


rule read_depth_filter:
    input:
        vcf_in = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"] + "/" + "{patient}.filtered.vcf.gz",
        tabix_in = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"] + "/" + "{patient}.filtered.vcf.gz.tbi"
    output:
        vcf_out = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"] + "/" + "{patient}.filtered_by_DP.vcf.gz"
    params:
        depth = config["params"]["read_depth_filter"]
    conda:
        "../envs/cyvcf2.yml"
    shell:
        """
        python3 scripts/filter_by_DP.py --input {input.vcf_in} --dp {params.depth}
        """

rule tabix_DP:
    input:
        vcf_in = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"] + "/" + "{patient}.filtered_by_DP.vcf.gz"
    output:
        tabix_out = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"] + "/" + "{patient}.filtered_by_DP.vcf.gz.tbi"
    conda:
        "../envs/tabix.yml"
    shell:
        """
        tabix -p vcf {input.vcf_in}
        """

rule select_calls:
    input:
        mutect2=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"] + "/" + "{patient}.filtered_by_DP.vcf.gz",
        strelka2=config["OUTPUT_FOLDER"] + config["OUTPUT_FOLDER"] + config["datadirs"]['VCF_out']+'/'+'{patient}_workflow'+'/results/variants/'+'variants.vcf.gz'
    output:
        overlap_vcf=config["OUTPUT_FOLDER"] + config["OUTPUT_FOLDER"] + config["datadirs"]['VCF_out']+'/'+'{patient}_overlap.vcf.gz'
    conda:
        "../envs/cyvcf2.yml"
    shell:
        """
        python3 scripts/overlap_calls.py --strelka {input.strelka2} --mutect {input.mutect2} --output {output.overlap_vcf}
        """

rule tabix_overlap:
    input:
        vcf_in = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_overlap.vcf.gz"
    output:
        tabix_out = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_overlap.vcf.gz.tbi"
    conda:
        "../envs/tabix.yml"
    shell:
        """
        tabix -p vcf {input.vcf_in}


