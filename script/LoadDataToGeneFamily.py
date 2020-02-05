#! /usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2019)

import os
import sys
import logging
from django.apps import apps
import datetime
from decimal import Decimal

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

def main(input_file, header):
    family_model = apps.get_model('expression_profiles', 'Gene_family')

    with open(input_file, 'r') as in_f:
        if header:
            next(in_f)
        for line in in_f:
            line = line.strip()
            tokens = line.split('\t')
            tmp = ['.'] * (9 - len(tokens))
            tokens.extend(tmp)
            data = family_model()
            update_dict = dict()
            human_gene_infos = models.Gene_info.objects.filter(species='Homo sapiens', status='activated')
            mouse_gene_infos = models.Gene_info.objects.filter(species='Mus musculus', status='activated')
            zebrafish_gene_infos = models.Gene_info.objects.filter(species='Danio rerio', status='activated')
            # human
            for gene_info in human_gene_infos:
                field_name = gene_info.name.lower() + '_id'
                gene_model = apps.get_model('expression_profiles', gene_info.name)
                genes = gene_model.objects.filter(ensembl_gene_id=tokens[1])
                if len(genes) > 0:
                    update_dict[field_name] = genes[0].id
            # mouse
            for gene_info in mouse_gene_infos:
                field_name = gene_info.name.lower() + '_id'
                gene_model = apps.get_model('expression_profiles', gene_info.name)
                genes = gene_model.objects.filter(ensembl_gene_id=tokens[3])
                if len(genes) > 0:
                    update_dict[field_name] = genes[0].id
            update_dict['mouse_support'] = tokens[5]
            # zebrafish
            for gene_info in zebrafish_gene_infos:
                field_name = gene_info.name.lower() + '_id'
                gene_model = apps.get_model('expression_profiles', gene_info.name)
                genes = gene_model.objects.filter(ensembl_gene_id=tokens[6])
                if len(genes) > 0:
                    update_dict[field_name] = genes[0].id
            update_dict['zebrafish_support'] = tokens[8]
            data.__dict__.update(update_dict)
            data.save()
if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Usage:
    %(prog)s -i merged_hsa_mmu_dps.tsv

    Input format(.tsv):
        1. family_ID
        2. human_ensembl_gene
        3. human_symbol
        4. mouse_ensembl_gene
        5. mouse_symbol
        6. mouse_support
        7. zebrafish_ensembl_gene
        8. zebrafish_symbol
        9. zebrafish_support
    """))

    parser.add_argument('-i', '--input', help='Input tsv file', required=True)
    parser.add_argument('-header', '--header', action='store_true', default=False, help='Ignore first line of the file')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    main(input_file=args.input, header=args.header)
