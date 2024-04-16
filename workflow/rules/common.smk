import pandas as pd
import os
import glob
from snakemake.utils import min_version

min_version("5.9.1")


configfile: "config/config_main.yaml"


# Load patient info.
# Note that this dataframe is accessed every time to determine the wildcards used
# at each step of the analysis while needed.

configpath = "config/config_main.yaml"

patients = pd.read_csv("patients.csv")["patient"]
units = pd.read_csv("units.csv").set_index(["patient"], drop=False)
units = units.sort_index()

slurm_logdir = config["slurm_log_dir"]
logpath = Path(slurm_logdir)
logpath.mkdir(parents=True, exist_ok=True)

bam_final_path = config["datadirs"]["BQSR"]
ref_fasta = config["resources"]["genome"]
ref_dict = ref_fasta.replace(".fa", ".dict")
HLA_path = config["OUTPUT_FOLDER"] + config["datadirs"]["HLA_typing"]
intervals_path = os.path.join(
    config["OUTPUT_FOLDER"] + config["datadirs"]["utils"], "interval-files"
)
vcfs_somatic_path = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF"]
vcfs_germline_path = config["OUTPUT_FOLDER"] + config["datadirs"]["VCF_germ"]

num_workers = 20


READ = ["1", "2"]


wildcard_constraints:
    patient="|".join(patients),


def get_fastq(wildcards):
    """Return a dict where keys are read1/2 and values are list of files"""
    return {
        "r1": units.loc[wildcards.patient, "fq1"],
        "r2": units.loc[wildcards.patient, "fq2"],
    }


def memory_for_gatk(gatk_mem: int):
    """Quick workaround to return string parameter for gatk"""
    as_str = f'--java-options "-Xmx{gatk_mem}g"'
    return as_str


def get_bam(wildcards):
    """
    Return list of bam files for rule bam_readcount
    """
    return [f"{bam_final_path}/{wildcards.patient}_recal.bam"]


def get_intervals():
    ints = []
    for i in range(num_workers):
        num_zeros = 4 - len(str(i))
        interval = "0" * num_zeros + str(i)
        ints.append(interval)
    return ints


def get_interval_files():
    ints = get_intervals()
    files = [i + "-scattered.interval_list" for i in ints]
    files = [os.path.join(intervals_path, f) for f in files]
    return files


def get_orientationbias_input(wildcards):
    intervals = get_intervals()
    files = [
        f"{vcfs_somatic_path}/{wildcards.patient}.{i}.f1r2.tar.gz" for i in intervals
    ]
    return files


def get_mergevcfs_input(wildcards):
    intervals = get_intervals()
    files = [
        f"{vcfs_somatic_path}/{wildcards.patient}.{i}.unfiltered.vcf.gz"
        for i in intervals
    ]
    return files


def get_mergestats_input(wildcards):
    intervals = get_intervals()
    files = [
        f"{vcfs_somatic_path}/{wildcards.patient}.{i}.unfiltered.vcf.gz.stats"
        for i in intervals
    ]
    return files


def read_hla(wildcards):
    """
    Read list HLA from allele_input_pvacseq.csv used as input for pvacseq
    """
    # the try/except is just a workaround as Snakemake evaluates functions before running rule, so this file
    # didn't exist before the rule HLA_typing is completed.
    hla_genotype = ""
    try:
        with open(
            f"{HLA_path}/{wildcards.patient}_output/allele_input_pvacseq.csv", "r"
        ) as hla_file:
            hla_genotype = hla_file.read()
    except FileNotFoundError:
        hla_genotype = "HLA-A*11:303"
    return hla_genotype


def sample_from_patient(df, patient_list, condition):
    samples = []
    for x in patient_list:
        samples.append(
            df.loc[(df["phenotype"] == condition) & (df["subject_id"] == x)].values[0]
        )
    return samples


interval_files = get_interval_files()
