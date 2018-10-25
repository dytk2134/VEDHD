#! /usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2018)

import os
import re
import sys
import subprocess
import logging
import vcf
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0.0'

def main(input, prefix):
    vcf_reader = vcf.Reader(open(input, 'r'))
    output_freq = prefix + '_freq'
    output_count = prefix + '_count'
    freq = open(output_freq, 'w')
    count = open(output_count, 'w')
    freq_header = ['chr', 'pos', 'ref' , 'alt','ExAC_AFR_Ref_freq', 'ExAC_AFR_Alt_freq','ExAC_AMR_Ref_freq', 'ExAC_AMR_Alt_freq','ExAC_EAS_Ref_freq', 'ExAC_EAS_Alt_freq','ExAC_FIN_Ref_freq', 'ExAC_FIN_Alt_freq','ExAC_NFE_Ref_freq', 'ExAC_NFE_Alt_freq','ExAC_OTH_Ref_freq', 'ExAC_OTH_Alt_freq','ExAC_SAS_Ref_freq', 'ExAC_SAS_Alt_freq' ]
    count_header = ['chr', 'pos', 'ref' , 'alt','ExAC_AFR_Ref_count', 'ExAC_AFR_Alt_count','ExAC_AMR_Ref_count', 'ExAC_AMR_Alt_count','ExAC_EAS_Ref_count', 'ExAC_EAS_Alt_count','ExAC_FIN_Ref_count', 'ExAC_FIN_Alt_count','ExAC_NFE_Ref_count', 'ExAC_NFE_Alt_count','ExAC_OTH_Ref_count', 'ExAC_OTH_Alt_count','ExAC_SAS_Ref_count', 'ExAC_SAS_Alt_count' ]
    freq.write('\t'.join(freq_header) + '\n')
    count.write('\t'.join(count_header) + '\n')
    population = ['AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
    for record in vcf_reader:
        for index, alt in enumerate(record.ALT):
            freq_outline = [record.CHROM, str(record.POS), str(record.REF), str(alt)]
            count_outline = [record.CHROM, str(record.POS), str(record.REF), str(alt)]
            #count_outline = [record.CHROM, str(record.POS), str(record.REF), str(alt), record.INFO['AN_AFR'][index], record.INFO['AC_AFR'], record.INFO['AN_AMR'][index], record.INFO['AC_AMR'],record.INFO['AN_EAS'][index], record.INFO['AC_EAS'], record.INFO['AN_FIN'][index], record.INFO['AC_FIN'], record.INFO['AN_NFE'][index], record.INFO['AC_NFE'], record.INFO['AN_OTH'][index], record.INFO['AC_OTH'],record.INFO['AN_SAS'][index], record.INFO['AC_SAS']]
            for pop in population:
                AN = 'AN_' + pop
                AC = 'AC_' + pop
                count_ref = record.INFO[AN] - record.INFO[AC][index]
                count_alt = record.INFO[AC][index]
                if record.INFO[AN] != 0:
                    freq_ref = count_ref/record.INFO[AN]
                    freq_alt = count_alt/record.INFO[AN]
                else:
                    freq_ref = 0
                    freq_alt = 0
                count_outline.extend([count_ref, count_alt])
                freq_outline.extend([freq_ref, freq_alt])
            freq.write('\t'.join(str(x) for x in freq_outline) + '\n')
            count.write('\t'.join(str(x) for x in count_outline) + '\n')
    freq.close()
    count.close()

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Usage:
    %(prog)s -input in.vcf -output out.tsv
    """))

    parser.add_argument('-in', '--input', type=str, help='Input vcf file', required=True)
    parser.add_argument('-prefix', '--prefix', type=str, help='Output tsv file', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    main(args.input, args.prefix)
