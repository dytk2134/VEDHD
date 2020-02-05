from django.db import models
from .models import GRCh37_Alleles, GRCh38_Alleles

# Create your models here.

class GRCh37_ESP_6500(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=256)
    Alt = models.CharField(max_length=256)
    EA_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EA_Ref_frep')
    EA_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EA_Alt_frep')
    AA_Ref_count = models.PositiveIntegerField(db_column='AA_Ref_count')
    AA_Alt_count = models.PositiveIntegerField(db_column='AA_Alt_count')

class GRCh38_ESP_6500(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=256)
    Alt = models.CharField(max_length=256)
    EA_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EA_Ref_frep')
    EA_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='EA_Alt_frep')
    AA_Ref_count = models.PositiveIntegerField(db_column='AA_Ref_count')
    AA_Alt_count = models.PositiveIntegerField(db_column='AA_Alt_count')

class GRCh37_ESP_variants(models.Model):
    allele_id = models.ForeignKey(GRCh37_Alleles, on_delete=models.CASCADE)
    grch37_esp_6500_id = models.ForeignKey(GRCh37_ESP_6500, on_delete=models.CASCADE)

class GRCh38_ESP_variants(models.Model):
    allele_id = models.ForeignKey(GRCh38_Alleles, on_delete=models.CASCADE)
    grch38_esp_6500_id = models.ForeignKey(GRCh38_ESP_6500, on_delete=models.CASCADE)
