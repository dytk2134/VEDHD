from django.db import models

# Create your models here.
class News(models.Model):
    release_date = models.DateField(blank=True, null=True)
    memo = models.TextField(blank=True, null=True)

class User_info(models.Model):
    user_id = models.CharField(max_length=32)
    request_date = models.DateField(blank=True, null=True)
    email = models.TextField(blank=True, null=True)
