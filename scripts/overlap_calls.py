#!/usr/bin/python3

from cyvcf2 import VCF
from cyvcf2 import Writer
import argparse

def overlap_calls(strelka_vcf: str, mutect_vcf: str, output_VCF: str):
    """
    Extract overlap calls between two VCF with a minimum depth required.
    """
    mutect_vcf = VCF(mutect_vcf)
    strelka_vcf = VCF(strelka_vcf)
    w = Writer(output_VCF, mutect_vcf, mode="wz")
    allowed_chrs = [str(x+1) for x in range(23)] + ['X','Y']    
    for variant in mutect_vcf:
        if str(variant.CHROM) in allowed_chrs:
            pos_as_str = f"{variant.CHROM}:{variant.POS}-{variant.POS + 1}"
            if len(list(strelka_vcf(pos_as_str))) > 0:
                for mut in strelka_vcf(pos_as_str):
                    # lower depth for strelka
                    if mut.FILTER == None and mut.gt_depths[0] > 20:
                        # check allele
                        if mut.ALT[0] == variant.ALT[0]:
                            w.write_record(variant)
                            break
            else:
                continue
                

def main():
    parser = argparse.ArgumentParser(description="Overlap VCFs")
    parser.add_argument("--strelka", help="Strelka VCF file")
    parser.add_argument("--mutect", help="Mutect VCF file")
    parser.add_argument("--output", help='Output VCF file')
    args=parser.parse_args()
    overlap_calls(args.strelka, args.mutect, args.output)

if __name__ == '__main__':
    main()
