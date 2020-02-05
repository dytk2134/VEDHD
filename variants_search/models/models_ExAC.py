from django.db import models
from .models import GRCh37_Alleles, GRCh38_Alleles

# Create your models here.

class GRCh37_ExAC_r1(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=400)
    Alt = models.CharField(max_length=450)
    AFR_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AFR_Ref_frep')
    AFR_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AFR_Alt_frep')
    AMR_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AMR_Ref_frep')
    AMR_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='AMR_Alt_frep')
    EAS_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EAS_Ref_frep')
    EAS_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EAS_Alt_frep')
    FIN_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='FIN_Ref_frep')
    FIN_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='FIN_Alt_frep')
    NFE_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='NFE_Ref_frep')
    NFE_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='NFE_Alt_frep')
    OTH_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='OTH_Ref_frep')
    OTH_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='OTH_Alt_frep')
    SAS_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='SAS_Ref_frep')
    SAS_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='SAS_Alt_frep')
    AFR_Ref_count = models.PositiveIntegerField(db_column='AFR_Ref_count')
    AFR_Alt_count = models.PositiveIntegerField(db_column='AFR_Alt_count')
    AMR_Ref_count = models.PositiveIntegerField(db_column='AMR_Ref_count')
    AMR_Alt_count = models.PositiveIntegerField(db_column='AMR_Alt_count')
    EAS_Ref_count = models.PositiveIntegerField(db_column='EAS_Ref_count')
    EAS_Alt_count = models.PositiveIntegerField(db_column='EAS_Alt_count')
    FIN_Ref_count = models.PositiveIntegerField(db_column='FIN_Ref_count')
    FIN_Alt_count = models.PositiveIntegerField(db_column='FIN_Alt_count')
    NFE_Ref_count = models.PositiveIntegerField(db_column='NFE_Ref_count')
    NFE_Alt_count = models.PositiveIntegerField(db_column='NFE_Alt_count')
    OTH_Ref_count = models.PositiveIntegerField(db_column='OTH_Ref_count')
    OTH_Alt_count = models.PositiveIntegerField(db_column='OTH_Alt_count')
    SAS_Ref_count = models.PositiveIntegerField(db_column='SAS_Ref_count')
    SAS_Alt_count = models.PositiveIntegerField(db_column='SAS_Alt_count')

class GRCh37_ExAC_variants(models.Model):
    allele_id = models.ForeignKey(GRCh37_Alleles, on_delete=models.CASCADE)
    grch37_exac_r1_id = models.ForeignKey(GRCh37_ExAC_r1, on_delete=models.CASCADE)