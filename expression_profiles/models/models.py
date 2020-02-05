from django.db import models
from .models_human import GRCh37_Genes, GRCh38_Genes
from .models_mouse import GRCm38_Genes
from .models_zebrafish import GRCz11_Genes
# Create your models here.

class Table_info(models.Model):
    name = models.CharField(max_length=64)
    species = models.CharField(max_length=64)
    version = models.CharField(max_length=32)
    tissue_id = models.TextField(blank=True, null=True)
    create_date = models.DateField(blank=True, null=True)
    modified_date = models.DateField(blank=True, null=True)
    status = models.CharField(max_length=32)
    memo = models.TextField(blank=True, null=True)

class Gene_info(models.Model):
    name = models.CharField(max_length=64)
    species = models.CharField(max_length=64)
    version = models.CharField(max_length=32)
    create_date = models.DateField(blank=True, null=True)
    modified_date = models.DateField(blank=True, null=True)
    status = models.CharField(max_length=32)
    memo = models.TextField(blank=True, null=True)

class Tissue_info(models.Model):
    tissue = models.CharField(max_length=64)
    abbreviation = models.CharField(max_length=64)

class Gene_family(models.Model):
    grch37_genes = models.ForeignKey(GRCh37_Genes, on_delete=models.CASCADE, null=True)
    grch38_genes = models.ForeignKey(GRCh38_Genes, on_delete=models.CASCADE, null=True)
    grcm38_genes = models.ForeignKey(GRCm38_Genes, on_delete=models.CASCADE, null=True)
    grcz11_genes = models.ForeignKey(GRCz11_Genes, on_delete=models.CASCADE, null=True)
    Mus_musculus_support = models.TextField(blank=True, null=True)
    Danio_rerio_support = models.TextField(blank=True, null=True)
