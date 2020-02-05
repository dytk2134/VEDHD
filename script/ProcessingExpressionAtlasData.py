#! /usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2019)

import os
import sys
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0.0'

def main(expression_info, output):
    header = list()
    with open(output, 'w') as out_f:
        with open(expression_info, 'r') as in_f:
            for line in in_f:
                line = line.strip()
                if line[0] == '#':
                    continue
                else:
                    tokens = line.split('\t')
                    if not tokens[0][0] == 'E':
                        header = list(tokens)
                    else:
                        gene_id = tokens[0]
                        tmp = [str()] * (len(header) - len(tokens))
                        tokens.extend(tmp)
                        for idx, value in enumerate(header[2:]):
                            tpm_value = tokens[idx + 2]
                            if tpm_value:
                                outline = [
                                    gene_id,
                                    value,
                                    tpm_value
                                ]
                                out_f.write('\t'.join(outline) + '\n')


if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    This script is for processing the data download from Expression Atlas.

    Input format:
    # Expression Atlas
    Gene ID	            Gene Name	brain	heart	kidney
    ENSMUSG00000000001	Gnai3	    44	    22	    74

    Output format:
    1. Ensembl Gene ID
    2. Tissue
    3. TPM Value

    Usage:
    %(prog)s -e E-GEOD-74747-query-results.tpms -o EGEOD_74747.txt
    """))

    parser.add_argument('-e', '--expression_info', type=str, help='Expression Data download from Expression Atlas', required=True)
    parser.add_argument('-o', '--output', type=str, help='The output filename', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    main(expression_info=args.expression_info, output=args.output)
