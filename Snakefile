
import glob
import os
import pandas as pd
import sys
from collections import Counter

configfile:
    "config_main.yaml"


def keep_paired_samples(patient_list):
    counts = Counter(patient_list)
    # keep only patients with matched tumor/normal
    to_keep = [patient for patient in patient_list if counts[patient] > 1]
    # escamoutage to keep insertion order. Works only on Python > 3.6 !
    return list(dict.fromkeys(to_keep).keys())


def sample_from_patient(df, patient_list, condition):
    samples = []
    for x in patient_list:
        samples.append(df[(df.phenotype == condition) &
                          (df.subject_id == x)].Sample.values[0])
    return samples


include:
    "rules/common.smk"

rule targets:
    input:        
        #expand(config["datadirs"]["VCF_out"] + "/" + "{patient}_ref_transcript_mismatch_reporter_filtered.vcf.gz", patient=patients),
        #expand(config["datadirs"]["VCF_out"]+"/"+"{patient}_vcf_readcount_annotator_FINAL.vcf", patient=patients)
        expand(config["datadirs"]["VCF_out"] + "/" + "{patient}_ref_transcript_mismatch_reporter.vcf.gz.tbi", patient=patients)
        
include:
    "rules/alignment.smk"
include:
    "rules/quantification.smk"
include:
    "rules/bam_cleaning.smk"
include:
    "rules/base_recalibration.smk"
include:
    "rules/HLA_typing.smk"
include:
    "rules/strelka.smk"
include:
    "rules/snv_calling.smk"
include:
    "rules/vatools.smk"
include:
    "rules/annotate_variants.smk"
include:
    "rules/decompose.smk"
include:
    "rules/bam_readcount.smk"
include:
    "rules/vcf_readcount_annotator_snp.smk"
include:
    "rules/vcf_readcount_annotator_indel.smk"
include:
    "rules/vcf_expression_annotator.smk"
include:
    "rules/ref_transcript_mismatch_reporter.smk"
include:
    "rules/pvacseq.smk"
