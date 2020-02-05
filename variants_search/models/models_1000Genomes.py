from django.db import models
from .models import GRCh37_Alleles, GRCh38_Alleles

# Create your models here.

class GRCh37_1000Genomes_Phase3(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=256)
    Alt = models.CharField(max_length=662)
    EAS_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EAS_Ref_frep')
    EAS_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EAS_Alt_frep')
    AMR_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AMR_Ref_frep')
    AMR_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AMR_Alt_frep')
    AFR_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AFR_Ref_frep')
    AFR_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AFR_Alt_frep')
    EUR_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EUR_Ref_frep')
    EUR_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EUR_Alt_frep')
    SAS_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='SAS_Ref_frep')
    SAS_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='SAS_Alt_frep')
    EAS_Ref_count = models.PositiveIntegerField(db_column='EAS_Ref_count')
    EAS_Alt_count = models.PositiveIntegerField(db_column='EAS_Alt_count')
    AMR_Ref_count = models.PositiveIntegerField(db_column='AMR_Ref_count')
    AMR_Alt_count = models.PositiveIntegerField(db_column='AMR_Alt_count')
    AFR_Ref_count = models.PositiveIntegerField(db_column='AFR_Ref_count')
    AFR_Alt_count = models.PositiveIntegerField(db_column='AFR_Alt_count')
    EUR_Ref_count = models.PositiveIntegerField(db_column='EUR_Ref_count')
    EUR_Alt_count = models.PositiveIntegerField(db_column='EUR_Alt_count')
    SAS_Ref_count = models.PositiveIntegerField(db_column='SAS_Ref_count')
    SAS_Alt_count = models.PositiveIntegerField(db_column='SAS_Alt_count')

class GRCh38_1000Genomes_Phase3(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=256)
    Alt = models.CharField(max_length=662)
    EAS_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EAS_Ref_frep')
    EAS_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EAS_Alt_frep')
    AMR_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AMR_Ref_frep')
    AMR_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AMR_Alt_frep')
    AFR_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AFR_Ref_frep')
    AFR_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AFR_Alt_frep')
    EUR_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EUR_Ref_frep')
    EUR_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EUR_Alt_frep')
    SAS_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='SAS_Ref_frep')
    SAS_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='SAS_Alt_frep')
    EAS_Ref_count = models.PositiveIntegerField(db_column='EAS_Ref_count')
    EAS_Alt_count = models.PositiveIntegerField(db_column='EAS_Alt_count')
    AMR_Ref_count = models.PositiveIntegerField(db_column='AMR_Ref_count')
    AMR_Alt_count = models.PositiveIntegerField(db_column='AMR_Alt_count')
    AFR_Ref_count = models.PositiveIntegerField(db_column='AFR_Ref_count')
    AFR_Alt_count = models.PositiveIntegerField(db_column='AFR_Alt_count')
    EUR_Ref_count = models.PositiveIntegerField(db_column='EUR_Ref_count')
    EUR_Alt_count = models.PositiveIntegerField(db_column='EUR_Alt_count')
    SAS_Ref_count = models.PositiveIntegerField(db_column='SAS_Ref_count')
    SAS_Alt_count = models.PositiveIntegerField(db_column='SAS_Alt_count')

class GRCh37_1000Genomes_variants(models.Model):
    allele_id = models.ForeignKey(GRCh37_Alleles, on_delete=models.CASCADE)
    grch37_1000genomes_phase3_id = models.ForeignKey(GRCh37_1000Genomes_Phase3, on_delete=models.CASCADE)

class GRCh38_1000Genomes_variants(models.Model):
    allele_id = models.ForeignKey(GRCh38_Alleles, on_delete=models.CASCADE)
    grch38_1000genomes_phase3_id = models.ForeignKey(GRCh38_1000Genomes_Phase3, on_delete=models.CASCADE)