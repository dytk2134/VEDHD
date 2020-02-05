from django.db import models
from .models import GRCh37_Alleles, GRCh38_Alleles

# Create your models here.

class GRCh37_TWB_NGS(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=256)
    Alt = models.CharField(max_length=256)
    TWN_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='TWN_Ref_frep')
    TWN_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='TWN_Alt_frep')
    TWN_Ref_count = models.PositiveIntegerField(db_column='TWN_Ref_count')
    TWN_Alt_count = models.PositiveIntegerField(db_column='TWN_Alt_count')

class GRCh37_TWB_GWG(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=256)
    Alt = models.CharField(max_length=256)
    TWN_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='TWN_Ref_frep')
    TWN_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='TWN_Alt_frep')
    TWN_Ref_count = models.PositiveIntegerField(db_column='TWN_Ref_count')
    TWN_Alt_count = models.PositiveIntegerField(db_column='TWN_Alt_count')

class GRCh37_TWB_variants(models.Model):
    allele_id = models.ForeignKey(GRCh37_Alleles, on_delete=models.CASCADE)
    grch37_twb_ngs_id = models.ForeignKey(GRCh37_TWB_NGS, on_delete=models.CASCADE)
    grch37_twb_gwg_id = models.ForeignKey(GRCh37_TWB_GWG, on_delete=models.CASCADE)