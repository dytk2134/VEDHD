#! /usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2019)

import os
import sys
import logging
import itertools

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0.0'

def main(gene_info, ensembl_biomark, output):
    ncbi_gene_info = {
        'GeneID': dict(),
        'Symbol': dict(),
        'Synonyms': dict(),
        'Ensembl_id': dict()
    }
    with open(gene_info, 'r') as in_f:
        next(in_f)
        for line in in_f:
            line = line.strip()
            tokens = line.split('\t')
            ncbi_gene_info['GeneID'][tokens[1]] = tokens # GeneID
            ncbi_gene_info['Symbol'][tokens[2]] = tokens # Symbol
            ncbi_gene_info['Synonyms'][tokens[4]] = tokens # Symbol
            dbXrefs = tokens[5].split('|')
            for dbXref in dbXrefs:
                if 'Ensembl:' in dbXref:
                    Ensembl_id = dbXref.replace('Ensembl:', '')
                    ncbi_gene_info['Ensembl_id'][Ensembl_id] = tokens
    All = 0
    MatchSymbol = 0
    Matchaleases = 0
    MatchEntrezID = 0
    MatchEnsemblID = 0
    NotMatch = 0
    gene_dict = dict()

    with open(output, 'w') as out_f:
        with open(ensembl_biomark, 'r') as in_f:
            next(in_f)
            for line in in_f:
                All += 1
                line = line.strip()
                tokens = line.split('\t')
                if len(tokens) != 8:
                    tmps = ['.'] * (8 - len(tokens))
                    tokens.extend(tmps)
                if tokens[0] not in gene_dict:
                    gene_dict[tokens[0]] = list()

                for idx, value in enumerate(tokens):
                    value = value.strip()
                    if not value or value == '-':
                        tokens[idx] = '.'
                # strand
                if tokens[4] == '1':
                    tokens[4] = '+'
                else:
                    tokens[4] = '-'

                outline = [
                    tokens[0], # Gene stable ID
                    tokens[6], # Gene Name
                    '.',       # NCBI gene name
                    '.',       # aliase
                    tokens[1], # Chromosome
                    tokens[2], # gene_start
                    tokens[3], # gene_end
                    tokens[4], # strand
                    tokens[5]  # gene description
                ]
                if tokens[0] in ncbi_gene_info['Ensembl_id']:
                    MatchEnsemblID += 1
                    outline[2] = ncbi_gene_info['Ensembl_id'][tokens[0]][2]
                    outline[3] = ncbi_gene_info['Ensembl_id'][tokens[0]][4].replace('|', ',')
                elif tokens[7] in ncbi_gene_info['GeneID']:
                    MatchEntrezID += 1
                    outline[2] = ncbi_gene_info['GeneID'][tokens[7]][2]
                    outline[3] = ncbi_gene_info['GeneID'][tokens[7]][4].replace('|', ',')
                elif tokens[6] in ncbi_gene_info['Symbol']:
                    MatchSymbol += 1
                    outline[2] = ncbi_gene_info['Symbol'][tokens[6]][2]
                    outline[3] = ncbi_gene_info['Symbol'][tokens[6]][4].replace('|', ',')
                else:
                    for ali in ncbi_gene_info['Synonyms']:
                        Synonyms = ali.split('|')
                        is_match = False
                        if tokens[6] in Synonyms:
                            Matchaleases += 1
                            outline[2] = ncbi_gene_info['Synonyms'][ali][2]
                            outline[3] = ncbi_gene_info['Synonyms'][ali][4].replace('|', ',')
                            is_match = True
                            break
                    if not is_match:
                        NotMatch += 1
                if not outline[2].strip() or outline[2] == '-':
                    outline[2] = '.'
                if not outline[3].strip() or outline[3] == '-':
                    outline[3] = '.'
                gene_dict[tokens[0]].append(outline)
        # Filter out duplicate
        for ensembl_id in gene_dict:
            gene_dict[ensembl_id].sort()
            # remove duplicate
            gene_dict[ensembl_id] = list(gene_dict[ensembl_id] for gene_dict[ensembl_id],_ in itertools.groupby(gene_dict[ensembl_id]))
            out_list = list()
            if len(gene_dict[ensembl_id]) == 1:
                out_list = list(gene_dict[ensembl_id])
            else:
                # Check if Ensembl gene name is same with NCBI gene name
                is_same_idx = set()
                for idx, tokens in enumerate(gene_dict[ensembl_id]):
                    if tokens[2] == '.' and tokens[3] == '.':
                        continue
                    else:
                        out_list.append(tokens)
                    if tokens[1] == tokens[2]:
                        is_same_idx.add(idx)
                if len(is_same_idx) != 0:
                    out_list = list()
                    for idx in is_same_idx:
                        out_list.append(gene_dict[ensembl_id][idx])
            if len(out_list) != 1:
                print(ensembl_id)
            for outline in out_list:
                out_f.write('\t'.join(outline) + '\n')
    print("Total row : %d" % (All))
    print("Symbol match : %d" % (MatchSymbol))
    print("Aleases match : %d" % (Matchaleases))
    print("Entrez ID match : %d" % (MatchEntrezID))
    print("Ensembl ID match : %d" % (MatchEnsemblID))
    print("Not match : %d" % (NotMatch))

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    This script is for generating the gene info table for VariED.

    Input:
    1. NCBI GENE_INFO (ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO)
    2. Ensembl BioMark table (tsv format):
        - Gene stable ID
        - Chromosome/scaffold Name
        - Gene start (bp)
        - Gene end (bp)
        - Strand
        - Gene description
        - Gene Name
        - NCBI gene ID (or EntrezGene ID)

    Usage:
    %(prog)s -g Homo_sapiens.gene_info -e Homo_sapiens_GRCh38p12.txt -o GRCh38_Genes.tsv
    """))

    parser.add_argument('-g', '--gene_info', type=str, help='NCBI GENE_INFO', required=True)
    parser.add_argument('-e', '--ensembl_biomark', type=str, help='Ensembl BioMark table (tsv format)', required=True)
    parser.add_argument('-o', '--output', type=str, help='The output filename of the gene info table for VariED', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    main(gene_info=args.gene_info, ensembl_biomark=args.ensembl_biomark, output=args.output)
