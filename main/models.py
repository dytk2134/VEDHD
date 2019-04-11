from django.db import models

# Create your models here.
class News(models.Model):
    release_date = models.DateField(blank=True, null=True)
    memo = models.TextField(blank=True, null=True)