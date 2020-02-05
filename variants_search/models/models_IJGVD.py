from django.db import models
from .models import GRCh37_Alleles, GRCh38_Alleles

# Create your models here.

class GRCh37_IJGVD_1KJPN(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=2)
    Alt = models.CharField(max_length=2)
    JPN_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='JPN_Ref_frep')
    JPN_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='JPN_Alt_frep')
    JPN_Ref_count = models.PositiveIntegerField(db_column='JPN_Ref_count')
    JPN_Alt_count = models.PositiveIntegerField(db_column='JPN_Alt_count')


class GRCh37_IJGVD_2KJPN(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=2)
    Alt = models.CharField(max_length=2)
    JPN_Ref_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='JPN_Ref_frep')
    JPN_Alt_frep = models.DecimalField(max_digits=5, decimal_places=4, db_column='JPN_Alt_frep')
    JPN_Ref_count = models.PositiveIntegerField(db_column='JPN_Ref_count')
    JPN_Alt_count = models.PositiveIntegerField(db_column='JPN_Alt_count')

class GRCh37_IJGVD_variants(models.Model):
    allele_id = models.ForeignKey(GRCh37_Alleles, on_delete=models.CASCADE)
    grch37_ijgvd_1kjpn_id = models.ForeignKey(GRCh37_IJGVD_1KJPN, on_delete=models.CASCADE)
    grch37_ijgvd_2kjpn_id = models.ForeignKey(GRCh37_IJGVD_2KJPN, on_delete=models.CASCADE)
