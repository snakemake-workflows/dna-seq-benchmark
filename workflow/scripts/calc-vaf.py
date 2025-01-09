from cyvcf2 import VCF

#vcf = "/Users/famke/01-pm4onco/osf-download/pipeline-results-of-imgag-data/qbic/strelka/tumor_5perc_vs_normal_5perc.strelka.somatic_snvs_VEP.ann.vcf" #snakemake.input
#indel = "/Users/famke/01-pm4onco/osf-download/pipeline-results-of-imgag-data/qbic/strelka/tumor_5perc_vs_normal_5perc.strelka.somatic_indels_VEP.ann.vcf.gz"

#
# SNVS
#
                                        
def get_snv_allele_freq(vcf):
    for variant in VCF(vcf):
        refCounts = variant.format(variant.REF + "U")
        altCounts = variant.format(variant.ALT[0] + "U")

        # TODO: check which value is the correct one from the matrix (this leads to many zero VAF)
        tier1RefCounts = refCounts[0, 0]
        tier1AltCounts = altCounts[0, 0]

        vaf = tier1AltCounts / (tier1AltCounts + tier1RefCounts)

        print(vaf)



#
# INDELs
#

def get_indel_allele_freq(vcf):
    for variant in VCF(vcf):
        tier1RefCounts = variant.format("TAR")[0,0]
        tier1AltCounts = variant.format("TIR")[0,0]

        vaf = tier1AltCounts / (tier1AltCounts + tier1RefCounts)

        print(vaf)
