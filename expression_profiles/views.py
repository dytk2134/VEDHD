from django.shortcuts import render, redirect
from django.template.defaulttags import register
from django.apps import apps
from expression_profiles import models
from django.views.decorators.csrf import csrf_protect
from django.urls import reverse
import uuid
import sys
import os
import collections
import re
import datetime
from django.contrib import messages
from django.db.models import Avg
# Create your views here.

delimiters = '\t',';',"\r\n",',',' '
regexPattern = '|'.join(map(re.escape,delimiters))

def multiple_replace(string, rep_dict):
    pattern = re.compile('|'.join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

@csrf_protect
def expression_profiles(request):
    # error message list
    messages = list()
    # get expression profiles table
    table_info_dict = dict()
    table_infos = models.Table_info.objects.filter(status='activated')
    for info in table_infos:
        if info.species not in table_info_dict:
            table_info_dict[info.species] = dict()
        table_info_dict[info.species][info.name] = list()
        for tissue_id in info.tissue_id.split(','):
            tissue_info = models.Tissue_info.objects.get(id=tissue_id)
            table_info_dict[info.species][info.name].append([tissue_info.tissue, tissue_info.abbreviation])
    if request.method == 'POST':
        job_id = str(uuid.uuid1())
        # processing user input gene list
        regexPattern = '|'.join(map(re.escape, delimiters))
        search_text = multiple_replace(request.POST.get('search_list'), {'\"': '', '\'': '', '`': '', '%': ''})
        search_text = search_text.upper()
        search_list = set(filter(None, re.split(regexPattern, search_text)))
        search_type = request.POST.get('SearchType')

        # the annotation table user choose. Currently, only allow one species
        # set default ref
        request_ref = request.POST.get('ref')
        if not request_ref:
            request_ref = 'GRCh37_Genes'
        annotation_table = apps.get_model('expression_profiles', request_ref)
        gene_info = models.Gene_info.objects.get(name=request_ref)
        annotation_title = ['User_input']
        reference_title = [f.name for f in annotation_table._meta.fields[1:]]
        annotation_title.extend(reference_title)

        duplicate_idx = list()
        duplicate_result = list()

        table_infos = models.Table_info.objects.filter(status='activated')
        table_field_dict = {
            'reference': collections.OrderedDict(),
            'other': collections.OrderedDict()
        }
        # expression profiles
        orthologs_title_dict = dict()
        for table_info in table_infos:
            tissue_list = request.POST.getlist(table_info.name)
            species = table_info.species.replace(' ', '_')
            if not tissue_info:
                continue
            elif tissue_list and table_info.species != gene_info.species:
                orthologs_title_dict[species] = ['ensembl_gene_name', 'ensembl_gene_id', species + '_support']

            field_list = list()
            for tissue in tissue_list:
                field_list.extend([tissue + '_TPM', tissue + '_Rank'])
                if table_info.species == gene_info.species:
                    # reference species
                    annotation_title.extend(['_'.join([table_info.name, species, tissue, 'TPM']), '_'.join([table_info.name, species, tissue, 'TPM'])])
                else:
                    # other species
                    orthologs_title_dict[species].extend(['_'.join([table_info.name, species, tissue, 'TPM']), '_'.join([table_info.name, species, tissue, 'TPM'])])
            if table_info.species == gene_info.species:
                table_field_dict['reference'][table_info.name] = list(field_list)
            else:
                table_field_dict['other'][table_info.name] = list(field_list)

        # get result of the target gene
        with open()
            line_number = 0
            for gene in search_list:
                gene_result = []
                if search_type == 'Symbol':
                    # ensembl
                    ensembl_result = annotation_table.filter(ensembl_gene_name=gene)
                    if ensembl_result:
                        gene_result = [ensembl_result[0]]
                    else:
                        # ncbi
                        ncbi_result = annotation_table.filter(NCBI_gene_name=gene)
                        if ncbi_result:
                            gene_result = [ncbi_result[0]]
                        else:
                            # aliases
                            aliases_result = annotation_table.objects.all().extra(where=['FIND_IN_SET("%s", aliases)' % (gene)])
                            if aliases_result:
                                gene_result = aliases_result
                                duplicate_list_id.append(aliases_result.values_list('id', flat=True))
                else:
                    ensembl_result = annotation_table.filter(ensembl_gene_id=gene)
                    if ensembl_result:
                        gene_result = [ensembl_result[0]]

                annotation_result = list()
                for gresult in gene_result:
                    result = [gene]
                    result.extend()
                    for ref_ex_table in table_field_dict['reference']:
                        expression_profile_table = apps.get_model('expression_profiles', ref_ex_table)
                        gene_id_field_name = request_ref.lower() + '_id__in'
                        for 
                        results = expression_profile_table.objects.filter(gene_id_field_name=alleles_id).values_list(*search_info_dict[search_info]['select_column'])
                if not annotation_table:
                    result = []
                    annotation_result.append(result)
                line_number += 1






        # search annotation table



        print(request.POST.getlist('ProteinAtlas_v18_1'))
    return render(request, 'expression_profiles/expression_profiles.html', locals())