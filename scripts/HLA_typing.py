#!/usr/bin/python3 

import sys

with open(sys.argv[1], 'r') as hla_input:
    HLAs = []
    for line in hla_input:
        if line.split('\t')[0] in ['HLA-A', 'HLA-B', 'HLA-C']:
            # check if homo/hetero zygous
            if line.split('\t')[1] == 2:
                HLAs.append(line.split('\t')[2].split(':')[:-1])
                HLAs.append(line.split('\t')[5].split(':')[:-1])
            else:
                HLAs.append(line.split('\t')[2].split(':')[:-1])
print(','.join(HLAs))
