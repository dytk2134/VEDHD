from django.db import models
from .models import GRCh37_Alleles, GRCh38_Alleles

# Create your models here.

class GRCh37_dbSNP_b151(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=1085)
    Alt = models.CharField(max_length=4291)
    rsID = models.CharField(max_length=32)

class GRCh37_dbSNP_b152(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=1085)
    Alt = models.CharField(max_length=4291)
    rsID = models.CharField(max_length=32)

class GRCh37_dbSNP_variants(models.Model):
    allele_id = models.ForeignKey(GRCh37_Alleles, on_delete=models.CASCADE)
    grch37_dbsnp_b151_id = models.ForeignKey(GRCh37_dbSNP_b151, on_delete=models.CASCADE)
    grch37_dbsnp_b152_id = models.ForeignKey(GRCh37_dbSNP_b152, on_delete=models.CASCADE)
