from django.db import models
from .models import GRCh37_Alleles, GRCh38_Alleles

# Create your models here.

class GRCh37_ClinVar_20190305(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=56)
    Alt = models.CharField(max_length=155)
    AlleleID = models.PositiveIntegerField()
    ClinicalSignificance = models.CharField(max_length=256)
    PhenotypeIDS = models.TextField()
    PhenotypeList = models.TextField()

class GRCh37_ClinVar_variants(models.Model):
    allele_id = models.ForeignKey(GRCh37_Alleles, on_delete=models.CASCADE)
    grch37_clinvar_20190305_id = models.ForeignKey(GRCh37_ClinVar_20190305, on_delete=models.CASCADE)
