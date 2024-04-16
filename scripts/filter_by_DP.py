#!/usr/bin/python3

from cyvcf2 import VCF
from cyvcf2 import Writer
import argparse


def filter_by_DP(vcf_file: str, read_depth: int):
    """
    Given a VCF file, it iterates over its records and output a vcf file where read depth
    of records is >= read depth

    Parameters
    ----------
    read_depth: int
        Chosen read depth
    vcf_file: str
        Path to a bgzipped & tab-indexed VCF file

    Returns
    ----------
    Filtered VCF file with "filtered.vcf.gz" suffix
    """
    transcr_VCF = VCF(vcf_file)
    output_VCF = vcf_file.split(".vcf.gz")[0] + "_by_DP.vcf.gz"
    w = Writer(output_VCF, transcr_VCF, mode="wz")
    for variant in transcr_VCF:
        if variant.gt_depths[0] >= read_depth:
            w.write_record(variant)
        else:
            continue


def main():
    parser = argparse.ArgumentParser(description="Filter VCF by chosen read depth")
    parser.add_argument("--input", help="Input VCF file")
    parser.add_argument(
        "--dp",
        type=int,
        default=30,
        help="Threshold to use for read depth filtering. Default 30",
    )
    args = parser.parse_args()
    if args.input == None:
        print("No input provided. Exiting..")
        exit()
    filter_by_DP(args.input, args.dp)


if __name__ == "__main__":
    main()
