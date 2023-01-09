rule pvacseq:
    input:
        vcf=config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"] + "/" + "{patient}_ref_transcript_mismatch_reporter_filtered.vcf.gz"
        # phaseVCF = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_out"]+"/"+"{patient}_annotated_phased.vcf.gz"
    params:
        samp_name = "{patient}",
        #hla_alleles = open(config["OUTPUT_FOLDER"] + config["datadirs"]["HLA_typing"]+'/'+'{wildcards.patient}_output/allele_input_pvacseq.csv', 'r').readline() ,
        hla_genotype=read_hla,
        #hla_alleles = '$(cat ' + config["OUTPUT_FOLDER"] + config["datadirs"]["HLA_typing"] + '/' + '{patient}_output' + '/' + 'allele_input_pvacseq.csv' + ')',
        e1 = config["params"]["pvacseq"]["e1"],
        mhc_tools = config["params"]["pvacseq"]["mhc_tools"],
        cores = config["params"]["pvacseq"]["t"],
        out_dir = config["OUTPUT_FOLDER"] + config["datadirs"]["pvacseq_out"],
        iedb_dir = config["resources"]["iedb_dir"]  
    output:
        out = config["OUTPUT_FOLDER"] + config["datadirs"]["pvacseq_out"] + "/{patient}_MHC_Class_I/" + "{patient}.all_epitopes.tsv"
    conda:
        "../envs/pvacseq.yaml"
    log:
        config["OUTPUT_FOLDER"] + config["datadirs"]["logs"]["pvacseq"] + '/{patient}.log'
    shell:
        """
        pvacseq run {input.vcf} {params.samp_name} {params.hla_genotype} \
        {params.mhc_tools} {params.out_dir} -e1 {params.e1} \
        --net-chop-method cterm --netmhc-stab -t {params.cores} \
        --iedb-install-directory {params.iedb_dir}
        """
