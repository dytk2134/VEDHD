from django.db import models
from .models import GRCh37_Alleles, GRCh38_Alleles

# Create your models here.

class GRCh37_REVEL(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=2)
    Alt = models.CharField(max_length=2)
    score = models.DecimalField(max_digits=5, decimal_places=4)

class GRCh37_REVEL_variants(models.Model):
    allele_id = models.ForeignKey(GRCh37_Alleles, on_delete=models.CASCADE)
    grch37_revel_id = models.ForeignKey(GRCh37_REVEL, on_delete=models.CASCADE)
