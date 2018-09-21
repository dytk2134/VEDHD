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

def main(input, output):
    vcf_reader = vcf.Reader(open(input, 'r'))
    for record in vcf_reader:
        print(record)

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Usage:
    %(prog)s -input in.vcf -output out.tsv
    """))

    parser.add_argument('-in', '--input', type=str, help='Input vcf file', required=True)
    parser.add_argument('-out', '--output', type=str, help='Output tsv file', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    main(args.input, args.output)
