# HLA typing will be performed using archasHLA (https://doi.org/10.1093/bioinformatics/btz474).
# This tool works basically by extracting reads on the sample's sixth chromosome and performing
# HLA genotyping via allele-specific quantification.
# Basically it's a porting of HISAT-genotype (https://www.biorxiv.org/content/10.1101/266197v1) strategy
# for RNA-seq. Now, its accuracy may vary a lot with different tissues as it's dependant on HLA locus'
# expression: authors declare very high performances (>90% accuracy, even in more complex tissues) and
# the software is already been used in many different studies.
# Besides that, its performance needs to be benchmarked against WES-WGS HLA genotyping, as it may be
# suffering in HNSC-specific samples.
#

# Rewrite rule to use only tumor bam

rule HLA_extract:
    input:
        bam = config["datadirs"]["BQSR_2"]+"/"+"{patient}_recal.pass2.bam"
        #bam = config["datadirs"]["mapped_reads"] +'/' + '{patient}_recal.pass2.bam.bam'
    output:
        fastq_1 = temp(config["datadirs"]["HLA_typing"] + '/' +'{patient}_output' + '/{patient}_recal.pass2' + '.extracted.1.fq.gz'),
        fastq_2 = temp(config["datadirs"]["HLA_typing"] + '/' +'{patient}_output' + '/{patient}_recal.pass2' + '.extracted.2.fq.gz')
    params:
        outdir = directory(config["datadirs"]["HLA_typing"] + '/' + '{patient}_output'),
        temp_dir = config['TEMP_DIR']
    conda:
        '../envs/arcasHLA.yml'
    threads:
        config["params"]["threads"]["arcasHLA"]
    log:
        config["datadirs"]["logs"]["arcasHLA"] + '/{patient}.log'
    shell:
        """
        arcasHLA extract {input.bam} -o {params.outdir} -t {threads} -v --temp {params.temp_dir}
        """
# Genotyping must be done before merging every genotype info inside a single .tsv file and converted
# using the 2 group notation (next rule). To do so, Snakemake offers the touch() function, which creates
# a file only after the job is concluded.
# Doing so, every job requiring the touched file in input will be executed only after the file it's
# created. This will be implemented only after seeing how arcasHLA creates genotype files.

rule HLA_genotype:
    input:
        fastq_1 = config["datadirs"]["HLA_typing"] + '/' +'{patient}_output' + '/{patient}_recal.pass2' + '.extracted.1.fq.gz',
        fastq_2 = config["datadirs"]["HLA_typing"] + '/' +'{patient}_output' + '/{patient}_recal.pass2' + '.extracted.2.fq.gz'
    output:
        genotype = config["datadirs"]["HLA_typing"] +'/' + '{patient}_output'+ '/{patient}_recal' + '.genotype.json'
    params:
        outdir = config["datadirs"]["HLA_typing"] +'/' + '{patient}_output',
        temp_dir = config['TEMP_DIR']
    conda:
        '../envs/arcasHLA.yml'
    threads:
        config["params"]["threads"]["arcasHLA"]
    log:
        config["datadirs"]["logs"]["arcasHLA"] + '/{patient}.log'
    shell:
        """
        arcasHLA genotype {input.fastq_1} {input.fastq_2} -o {params.outdir} -t {threads} --temp {params.temp_dir} 
        """

rule HLA_merge:
    input:
        genotype = config["datadirs"]["HLA_typing"] +'/' + '{patient}_output'+ '/{patient}_recal' + '.genotype.json'
    output:
        genotype_out = config["datadirs"]["HLA_typing"] + '/' +'{patient}_output' + '/{patient}_recal' + '.genotypes.tsv'
    params:
        outdir = config["datadirs"]["HLA_typing"] +'/' + '{patient}_output',
        name_out = '/{patient}_recal'
    conda:
        '../envs/arcasHLA.yml'
    threads:
        config["params"]["threads"]["arcasHLA"]
    log:
        config["datadirs"]["logs"]["arcasHLA"] + '/{patient}.log'
    shell:
        """
        arcasHLA merge --run {params.name_out} --i {params.outdir} --o {params.outdir} -v
        """

rule HLA_convert_resolution:
    input:
        genotype = config["datadirs"]["HLA_typing"] + '/' +'{patient}_output' + '/{patient}_recal' + '.genotypes.tsv'
    output:
        genotype_converted = config["datadirs"]["HLA_typing"] + '/' + '{patient}_output' + '/' + 'resolution_genotypes.tsv'
    conda:
        '../envs/arcasHLA.yml'
    log:
        config["datadirs"]["logs"]["arcasHLA"] + '/{patient}.log'
    shell:
        """
        arcasHLA convert --resolution 2 {input.genotype} -o {output.genotype_converted}
        """


rule extract_hla:
    input: 
        config["datadirs"]["HLA_typing"] + '/' + '{patient}_output' + '/' + 'resolution_genotypes.tsv'
    output: 
        config["datadirs"]["HLA_typing"] + '/' + '{patient}_output' + '/' + 'allele_input_pvacseq.csv'
    params:
        patient = '{patient}'
    log:
        config["datadirs"]["logs"]["arcasHLA"] + '/{patient}.log'
    run:
        hla_types = extract_hla(params.patient)
