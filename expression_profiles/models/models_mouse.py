from django.db import models

# Create your models here.

class GRCm38_Genes(models.Model):
    ensembl_gene_id = models.CharField(max_length=20)
    ensembl_gene_name = models.CharField(max_length=32)
    NCBI_gene_name = models.CharField(max_length=32)
    aliases = models.CharField(max_length=300)
    chromosome = models.CharField(max_length=2)
    gene_start = models.PositiveIntegerField()
    gene_end = models.PositiveIntegerField()
    strand = models.CharField(max_length=1)
    gene_description = models.TextField(null=True)

    class Meta:
        indexes = [
            models.Index(fields=['ensembl_gene_id']),
            models.Index(fields=['chromosome', 'gene_start', 'gene_end']),
            models.Index(fields=['ensembl_gene_name']),
            models.Index(fields=['NCBI_gene_name']),
            models.Index(fields=['aliases']),
        ]

class EGEOD_74747(models.Model):
    grcm38_genes = models.ForeignKey(GRCm38_Genes, on_delete=models.CASCADE, null=True)
    Brain_Rank = models.CharField(max_length=32, null=True)
    Brain_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Heart_Rank = models.CharField(max_length=32, null=True)
    Heart_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Kidney_Rank = models.CharField(max_length=32, null=True)
    Kidney_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Liver_Rank = models.CharField(max_length=32, null=True)
    Liver_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Lung_Rank = models.CharField(max_length=32, null=True)
    Lung_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Spleen_Rank = models.CharField(max_length=32, null=True)
    Spleen_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Testis_Rank = models.CharField(max_length=32, null=True)
    Testis_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Thymus_Rank = models.CharField(max_length=32, null=True)
    Thymus_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Skin_Rank = models.CharField(max_length=32, null=True)
    Skin_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
