from django.db import models


# Create your models here.
class Table_info(models.Model):
    name = models.CharField(max_length=64)
    species = models.CharField(max_length=64)
    version = models.CharField(max_length=32)
    population_id = models.TextField(blank=True, null=True)
    create_date = models.DateField(blank=True, null=True)
    modified_date = models.DateField(blank=True, null=True)
    status = models.CharField(max_length=32)
    memo = models.TextField(blank=True, null=True)

class Population_info(models.Model):
    population = models.CharField(max_length=64)
    abbreviation = models.CharField(max_length=64)

class GRCh37_Alleles(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=1085)
    Alt = models.CharField(max_length=4291)

    class Meta:
        indexes = [
            models.Index(fields=['Chr', 'Pos'])
        ]

class GRCh38_Alleles(models.Model):
    Chr = models.CharField(max_length=2)
    Pos = models.PositiveIntegerField()
    Ref = models.CharField(max_length=1085)
    Alt = models.CharField(max_length=4291)

    class Meta:
        indexes = [
            models.Index(fields=['Chr', 'Pos'])
        ]
