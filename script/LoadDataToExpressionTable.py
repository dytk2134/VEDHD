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

def get_rank(num, num_set):
    num_list = list(num_set)
    num_list_sorted = sorted(num_list, reverse=True)
    for idx, value in enumerate(num_list_sorted):
        if num == value:
            return str(idx+1) + '/' + str(len(num_list))

def main(model, input, header, species, data_version, memo):
    # expression profiles
    expression_model = apps.get_model('expression_profiles', model)
    expression_dict = dict()
    tissue_dict = dict()
    tissue_ids = set()
    rank_dict = dict()
    with open(input, 'r') as in_f:
        if header:
            next(in_f)
        for line in in_f:
            line = line.strip()
            tokens = line.split('\t')
            if tokens[0] not in expression_dict:
                expression_dict[tokens[0]] = dict()
            tokens[1] = tokens[1].capitalize()
            tissue = tokens[1].replace(' ', '_')
            expression_dict[tokens[0]][tissue] = Decimal(tokens[2])
            if tissue not in rank_dict:
                rank_dict[tissue] = set()
                tissue_dict[tissue] = tokens[1]
                # check if this tissue in tissue info
                tissue_infos = models.Tissue_info.objects.filter(abbreviation=tissue)
                if len(tissue_infos) == 0:
                    tissue_info = models.Tissue_info(tissue=tokens[1], abbreviation=tissue)
                    tissue_info.save()
                else:
                    tissue_info = tissue_infos[0]
                tissue_ids.add(tissue_info.id)
            rank_dict[tissue].add(Decimal(tokens[2]))

    # current time
    now = datetime.datetime.now()
    # check if model in the table
    table_infos = models.Table_info.objects.filter(name=model)
    if len(table_infos) > 0:
        table_info = table_infos[0]
        table_info.modified_date = now
    else:
        table_info = models.Table_info()
        table_info.create_date = now
        table_info.modified_date = now
        table_info.name = model
    table_info.species = species
    if not data_version:
        data_version = '.'
    table_info.version = data_version
    table_info.status = 'activated'
    if not memo:
        memo = '.'
    table_info.memo = memo
    table_info.tissue_id = ','.join(map(str, list(tissue_ids)))
    table_info.save()

    # get gene_info table
    gene_infos = models.Gene_info.objects.filter(species=species, status='activated')
    if len(gene_infos) == 0:
        logger.error('Failed to find any activated Gene info table of %s' % (species))
        sys.exit(1)
    id_fields = list()
    gene_models = list()
    for gene_info in gene_infos:
        field_name = gene_info.name.lower()
        id_fields.append(field_name)
        gene_model = apps.get_model('expression_profiles', gene_info.name)
        gene_models.append(gene_model)

    # expression profiles table
    tissue_set = set()
    for gene_id in expression_dict:
        data = expression_model()
        all_fields = data._meta.fields
        update_dict = dict()
        ignore_field = set()
        for field in all_fields:
            # tissue
            if field.name == 'id':
                continue
            elif field.name not in id_fields:
                tmp = field.name.split('_')[:-1]
                tissue = '_'.join(tmp)
                tissue_set.add(tissue)
                if tissue in expression_dict[gene_id]:
                    if 'Rank' in field.name:
                        # Rank
                        update_dict[field.name] = get_rank(expression_dict[gene_id][tissue], rank_dict[tissue])
                    else:
                        # TPM
                        update_dict[field.name] = expression_dict[gene_id][tissue]
            else:
                ignore_field.add(field.name)
        not_found_fields = ignore_field - set(id_fields)
        if len(not_found_fields) != 0:
            logger.error('Failed to find fields %s' % (str(not_found_fields)))
            sys.exit(1)

        # gene_id
        for idx, gene_model in enumerate(gene_models):
            id_field = id_fields[idx]
            if id_field in ignore_field:
                genes = gene_model.objects.filter(ensembl_gene_id=gene_id)
                if len(genes) > 0:
                    update_dict[id_field] = genes[0].id

        data.__dict__.update(update_dict)
        data.save()
        sys.exit(1)
    # check tissue
    not_found_tissue = set(tissue_dict.keys()) - tissue_set
    if len(not_found_tissue) != 0:
        logger.warning('Failed to find tissue fields: %s' % (str(not_found_tissue)))

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Usage:
    %(prog)s -m ProteinAtlas_v18_1 -input ProteinAtlas_v18_1.txt -s \"Homo sapiens\" -ver \"v18.1\" -memo \"Download from The Protein Atlas\"

    Input format(.tsv):
        1. ensembl_gene_id
        2. tissue
        3. Value(TPM)
    """))

    parser.add_argument('-m', '--model', help='Model name of gene expression table', required=True)
    parser.add_argument('-i', '--input', help='Input tsv file', required=True)
    parser.add_argument('-header', '--header', action='store_true', default=False, help='Ignore first line of the file')
    parser.add_argument('-species', '--species', help='Species. Format: [Genus] [species]', required=True)
    parser.add_argument('-data_version', '--data_version', help='The version of the input data.')
    parser.add_argument('-memo', '--memo', help='memo')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    main(args.model, args.input, args.header, args.species, args.data_version, args.memo)
