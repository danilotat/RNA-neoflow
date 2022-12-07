import sys


class GTF_element(object):
    def __init__(self,chromosome,source,feature_type,start,end,score,strand,phase,attributes):
        self.chromosome=chromosome
        self.source=source
        self.feature_type=feature_type
        self.start=start 
        self.end=end 
        self.score=score
        self.strand=strand
        self.phase=phase
        self.attributes=attributes 

    def parse_attributes(self):
        ''' Parse the 9th column of GTF into a dictionary of attributes.'''
        f_dict = {}
        kvp = self.attributes.split(';')[:-1]
        for i in kvp:
            ri = i.replace(' "', '="')
            try:
                k, v = ri.split('=')
                rk = k.replace(' ', '')
                f_dict[rk] = v
            except ValueError:
                print('Malformed file.\nThis line isn\'t good. %s' % kvp)
                exit()
        return f_dict



def main(gtf):
    with open(gtf,'r') as gtf_file:
        for line in gtf_file:
            if line.startswith('#'):
                continue
            else:
                entry = GTF_element(*line.rstrip().split('\t'))
                if entry.feature_type == "exon":
                    entry_dict = entry.parse_attributes()
                    if entry_dict["gene_biotype"] == '"protein_coding"':
                        print(f"{entry.chromosome}\t{int(entry.start) - 1}\t{entry.end}")

if __name__ == "__main__":
    main(sys.argv[1])