from django.db import models

# Create your models here.

class GRCz11_Genes(models.Model):
    ensembl_gene_id = models.CharField(max_length=20)
    ensembl_gene_name = models.CharField(max_length=32)
    NCBI_gene_name = models.CharField(max_length=32)
    aliases = models.CharField(max_length=300)
    chromosome = models.CharField(max_length=2)
    gene_start = models.PositiveIntegerField()
    gene_end = models.PositiveIntegerField()
    strand = models.CharField(max_length=1)
    gene_description = models.TextField(blank=True, null=True)

    class Meta:
        indexes = [
            models.Index(fields=['ensembl_gene_id']),
            models.Index(fields=['chromosome', 'gene_start', 'gene_end']),
            models.Index(fields=['ensembl_gene_name']),
            models.Index(fields=['NCBI_gene_name']),
            models.Index(fields=['aliases']),
        ]

class PRJNA293022(models.Model):
    grcz11_genes = models.ForeignKey(GRCz11_Genes, on_delete=models.CASCADE, null=True)
    Atrium_Rank = models.CharField(max_length=32, null=True)
    Atrium_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Bulbus_arteriosus_Rank = models.CharField(max_length=32, null=True)
    Bulbus_arteriosus_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Ventricle_Rank = models.CharField(max_length=32, null=True)
    Ventricle_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)

