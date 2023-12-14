from django.db import models

# Create your models here.

class UserUploadedMatrix(models.Model):
  matrix=models.TextField()