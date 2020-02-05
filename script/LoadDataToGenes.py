#! /usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2019)

import os
import sys
import logging
from django.apps import apps
import datetime

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings.settings')

from django.conf import settings
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
import expression_profiles.models as models

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0.0'

def main(model_name, input_file, header, species, data_version, memo):
    gene_model = apps.get_model('expression_profiles', model_name)
    with open(input_file, 'r') as in_f:
        if header:
            next(in_f)
        for line in in_f:
            line = line.strip()
            tokens = line.split('\t')
            gene = gene_model()
            gene.ensembl_gene_id = tokens[0]
            gene.ensembl_gene_name = tokens[1]
            gene.NCBI_gene_name = tokens[2]
            gene.aliases = tokens[3]
            gene.chromosome = tokens[4]
            gene.gene_start = tokens[5]
            gene.gene_end = tokens[6]
            gene.strand = tokens[7]
            gene.gene_description = tokens[8]
            gene.save()

    # current time
    now = datetime.datetime.now()
    gene_infos = models.Gene_info.objects.filter(name=model_name)
    if len(gene_infos) > 0:
        gene_info = gene_infos[0]
    else:
        gene_info = models.Gene_info()
        gene_info.name = model_name
        gene_info.create_date = now
    gene_info.modified_date = now
    gene_info.species = species
    if not data_version:
        data_version = '.'
    gene_info.version = data_version
    gene_info.status = 'activated'
    if not memo:
        memo = '.'
    gene_info.memo = memo
    gene_info.save()
if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Usage:
    %(prog)s -m GRCh37_Genes -input GRCh37_Genes.txt

    Input format(.tsv):
        1. ensembl_gene_id
        2. ensembl_gene_name
        3. NCBI_gene_name
        4. aliases
        5. chromosome
        6. gene_start
        7. gene_end
        8. strand
        9. gene_description
    """))

    parser.add_argument('-m', '--model', help='Model name of [version]_Genes table', required=True)
    parser.add_argument('-i', '--input', help='Input tsv file', required=True)
    parser.add_argument('-header', '--header', action='store_true', default=False, help='Ignore first line of the file')
    parser.add_argument('-species', '--species', help='Species. Format: [Genus] [species]', required=True)
    parser.add_argument('-data_version', '--data_version', help='The version of the input data.')
    parser.add_argument('-memo', '--memo', help='memo')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    main(args.model, args.input, args.header, args.species, args.data_version, args.memo)
