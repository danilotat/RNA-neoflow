import cyvcf2
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys

########### PRIORS ###########
pis = np.array(["0.00011384","0.00011123","0.0001134","0.00011447","0.00011315","0.00011356",
    "0.00011308","0.00011137","0.0001118","0.0001147","0.00011516","0.0001114","0.00011502",
    "0.00011414","0.00011204","0.00011251","0.00011369","0.00011515","0.00011314","0.00011102",
    "0.000114","0.00011275","0.00011452","0.00011382","0.00011344","0.00011177","0.0001122",
    "0.00011184","0.00011313","0.00011059","0.00011443","0.00011234","0.0001115","0.00011236",
    "0.00011078","0.00011357","0.00011404","0.00011232","0.00011393","0.00011206","0.00011411",
    "0.00011236","0.0001111","0.00011471","0.0001135","0.00011332","0.00011662","0.0001138",
    "0.00011056","0.00011303","0.00011252","0.00011169","0.00011245","0.00010985","0.00011473",
    "0.00011102","0.00011303","0.00011332","0.00011057","0.00011218","0.00011225","0.00011305",
    "0.00011611","0.0001131","0.00011185","0.00011408","0.00011381","0.00011215","0.00011511",
    "0.000113"], dtype=float)

alpha_exonMut, beta_exonMut = (0.8651947104981104, 185.90089630396128)
alpha_MutRatio, beta_MutRatio = (0.47724493362038634, 20.01410347598762)
alpha_pi, beta_pi = (6570.723319753332, 58154364.249282956)

##################################


def get_avg_freq(variant: cyvcf2.Variant):
    infodict = dict(variant.INFO)
    ANs = []
    AFs = []
    for AF, AN in zip(
        [x for x in infodict.keys() if 'AF' in x],
        [x for x in infodict.keys() if 'AN' in x]):
        ANs.append(infodict[AN])
        AFs.append(infodict[AF])
    # drop nan in both arrays
    ANs = [x for x in ANs if isinstance(x, float) or isinstance(x, int)]
    AFs = [x for x in AFs if isinstance(x, float) or isinstance(x, int)]
    if len(ANs) == len(AFs) and len(AFs) > 0:
        try:
            return np.average(AFs, weights=ANs)
        except ZeroDivisionError:
            return 0
    else:
        return np.nan

def phred_to_prob(phred: float):
    return 10**(-phred/10)


def prob_to_phred(prob: float):
    return -10*np.log10(prob)

class variantProb():
    """
    This class extends the cyvcf2.Variant class to include the probability of being a germline variant.
    """
    def __init__(self, variant: cyvcf2.Variant, pis: np.array):
        self.variant = variant
        self.pis = pis
        self.fixed_pi = stats.beta.rvs(alpha_pi, beta_pi)
        self.avg_freq = self._get_avg_freq()
        self.germ_LTs = self._assign_germ_LT()
        self.germ_prob = self._germProb()

    def _get_avg_freq(self):
        infodict = dict(self.variant.INFO)
        ANs = []
        AFs = []
        for AF, AN in zip(
            [x for x in infodict.keys() if 'AF' in x],
            [x for x in infodict.keys() if 'AN' in x]):
            ANs.append(infodict[AN])
            AFs.append(infodict[AF])
        # drop nan in both arrays
        ANs = [x for x in ANs if isinstance(x, float) or isinstance(x, int)]
        AFs = [x for x in AFs if isinstance(x, float) or isinstance(x, int)]
        if len(ANs) == len(AFs) and len(AFs) > 0:
            try:
                return np.average(AFs, weights=ANs)
            except ZeroDivisionError:
                return 0
        else:
            # when we have no AF or AN, we use the priors as defined in the Mutect2 paper
            return stats.beta.rvs(0.01, 10)
        
        
        
    def _assign_germ_LT(self) -> list:
        """
        Goal of this function is to assign likelihood for the only three main cases.
        1. Ref homo
        2. Heterozygous
        3. Alt homo
        """
        LTs = []
        try:
            for gt_LT in (self.variant.gt_phred_ll_homref[0],self.variant.gt_phred_ll_het[0],self.variant.gt_phred_ll_homalt[0]):
                LTs.append(phred_to_prob(gt_LT))
            return LTs
        except IndexError:
            return [0,0,0]
        
    def __str__(self):
        return str(self.variant)


    def _germProb(self) -> tuple:
        p_S = self.fixed_pi
        p0_1 = phred_to_prob(self.germ_LTs[1]) * self.avg_freq * (1 - self.avg_freq) 
        p1_1 = self.avg_freq**2 * phred_to_prob(self.germ_LTs[-1]) 
        # get normalized probabilities
        post0_1 = prob_to_phred(p0_1 / (p_S + p0_1 + p1_1))
        post1_1 = prob_to_phred(p1_1 / (p_S + p0_1 + p1_1))
        p_S = prob_to_phred(p_S / (p_S + p0_1 + p1_1))
        return (post0_1, post1_1, p_S)
        

def writeGermProb(input_vcf: str, output_vcf: str):
    """
    This function writes the germProb to the VCF file
    """
    vcf_file = cyvcf2.VCF(input_vcf, threads=8)
    vcf_file.add_info_to_header(
        {'ID': 'hetProb', 'Description': 'PHRED score of germline for alt heterozygosity given population allele frequency', 'Type': 'Float', 'Number': "1"})
    vcf_file.add_info_to_header(
        {'ID': 'homaltProb', 'Description': 'PHRED score of germline for alt homozygosity given population allele frequency', 'Type': 'Float', 'Number': "1"})
    vcf_file.add_info_to_header(
        {'ID': 'somProb', 'Description': 'PHRED score of somatic event for all genotypes given population allele frequency', 'Type': 'Float', 'Number': "1"}
    )
    outfile = cyvcf2.Writer(output_vcf, tmpl=vcf_file, mode='wz')
    for variant in vcf_file:
        variant_adj = variantProb(variant, pis)
        variant.INFO['hetProb'] = variant_adj.germ_prob[0]
        variant.INFO['homaltProb'] = variant_adj.germ_prob[1]
        variant.INFO['somProb'] = variant_adj.germ_prob[2]
        outfile.write_record(variant)
    outfile.close()
    


if __name__ == "__main__":
    writeGermProb(sys.argv[1], sys.argv[2])