from django.db import models

# Create your models here.

class GRCh37_Genes(models.Model):
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

class GRCh38_Genes(models.Model):
    ensembl_gene_id = models.CharField(max_length=15)
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

class ProteinAtlas_v18_1(models.Model):
    grch37_genes = models.ForeignKey(GRCh37_Genes, on_delete=models.CASCADE, null=True)
    grch38_genes = models.ForeignKey(GRCh38_Genes, on_delete=models.CASCADE, null=True)
    Adipose_tissue_Rank = models.CharField(max_length=32, null=True)
    Adipose_tissue_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Adrenal_gland_Rank = models.CharField(max_length=32, null=True)
    Adrenal_gland_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Appendix_Rank = models.CharField(max_length=32, null=True)
    Appendix_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Bone_marrow_Rank = models.CharField(max_length=32, null=True)
    Bone_marrow_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Breast_Rank = models.CharField(max_length=32, null=True)
    Breast_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Cerebral_cortex_Rank = models.CharField(max_length=32)
    Cerebral_cortex_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Cervix_uterine_Rank = models.CharField(max_length=32, null=True)
    Cervix_uterine_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Colon_Rank = models.CharField(max_length=32, null=True)
    Colon_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Duodenum_Rank = models.CharField(max_length=32, null=True)
    Duodenum_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Endometrium_Rank = models.CharField(max_length=32, null=True)
    Endometrium_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Epididymis_Rank = models.CharField(max_length=32, null=True)
    Epididymis_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Esophagus_Rank = models.CharField(max_length=32, null=True)
    Esophagus_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Fallopian_tube_Rank = models.CharField(max_length=32, null=True)
    Fallopian_tube_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Gallbladder_Rank = models.CharField(max_length=32, null=True)
    Gallbladder_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Heart_muscle_Rank = models.CharField(max_length=32, null=True)
    Heart_muscle_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Kidney_Rank = models.CharField(max_length=32, null=True)
    Kidney_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Liver_Rank = models.CharField(max_length=32, null=True)
    Liver_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Lung_Rank = models.CharField(max_length=32, null=True)
    Lung_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Lymph_node_Rank = models.CharField(max_length=32, null=True)
    Lymph_node_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Ovary_Rank = models.CharField(max_length=32, null=True)
    Ovary_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Pancreas_Rank = models.CharField(max_length=32, null=True)
    Pancreas_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Parathyroid_gland_Rank = models.CharField(max_length=32, null=True)
    Parathyroid_gland_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Placenta_Rank = models.CharField(max_length=32, null=True)
    Placenta_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Prostate_Rank = models.CharField(max_length=32, null=True)
    Prostate_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Rectum_Rank = models.CharField(max_length=32, null=True)
    Rectum_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Salivary_gland_Rank = models.CharField(max_length=32, null=True)
    Salivary_gland_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Seminal_vesicle_Rank = models.CharField(max_length=32, null=True)
    Seminal_vesicle_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Skeletal_muscle_Rank = models.CharField(max_length=32, null=True)
    Skeletal_muscle_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Skin_Rank = models.CharField(max_length=32, null=True)
    Skin_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Small_intestine_Rank = models.CharField(max_length=32, null=True)
    Small_intestine_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Smooth_muscle_Rank = models.CharField(max_length=32, null=True)
    Smooth_muscle_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Spleen_Rank = models.CharField(max_length=32, null=True)
    Spleen_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Stomach_Rank = models.CharField(max_length=32, null=True)
    Stomach_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Testis_Rank = models.CharField(max_length=32, null=True)
    Testis_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Thyroid_gland_Rank = models.CharField(max_length=32, null=True)
    Thyroid_gland_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Tonsil_Rank = models.CharField(max_length=32, null=True)
    Tonsil_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
    Urinary_bladder_Rank = models.CharField(max_length=32, null=True)
    Urinary_bladder_TPM = models.DecimalField(max_digits=11, decimal_places=2, null=True)
