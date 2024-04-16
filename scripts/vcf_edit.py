#!/usr/bin/python3

import sys

test_set = [
    "chr1",
    "chr10_GL383545v1_alt",
    "chr10_GL383546v1_alt",
    "chrUn_KI270420v1",
    "chrUn_KI270422v1",
    "chr10_KI270824v1_alt",
    "chr10_KI270825v1_alt",
    "chr11",
    "chr11_GL383547v1_alt",
    "chr11_JH159136v1_alt",
]


def avail_chr(chr_list_file):
    chr_genome_name = []
    with open(chr_list_file, "r") as chr_list:
        for line in chr_list:
            chr_genome_name.append(line.rstrip())
    return chr_genome_name


def ncbi_to_ensembl(chr: str) -> str:
    if "HLA" in chr:
        return str(chr)
    else:
        if chr.count("_") == 0:
            if "chr" in chr:
                chr_conv = chr.replace("chr", "")
                return str(chr_conv)
            else:
                print("Invalid chromosome name!")
                sys.exit()
        elif 3 > chr.count("_") >= 1:
            chr = chr.split("_")[1].replace("v", ".")
            return str(chr)
        else:
            print("Invalid character, exiting..")
            sys.exit()


if __name__ == "__main__":
    good_res = avail_chr(sys.argv[1])
    with open(sys.argv[2], "r") as contig_list:
        for line in contig_list:
            ens_record = ncbi_to_ensembl(line.rstrip())
            if ens_record in good_res:
                print("{}\t{}".format(line.rstrip(), ens_record))
            else:
                print("{}\t{}".format(line.rstrip(), line.rstrip()))
