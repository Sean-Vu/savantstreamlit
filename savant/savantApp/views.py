from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt  # Add this import
from .models import UserUploadedMatrix
import json

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

@csrf_exempt  # Add this decorator to disable CSRF protection for simplicity in this example
def upload_matrix(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body.decode('utf-8'))
            matrix_content = data.get('matrix')
            print(matrix_content)
            
            user_uploaded_matrix = UserUploadedMatrix(matrix=matrix_content)
            user_uploaded_matrix.save()

            return JsonResponse({'status': 'success'})
        except json.JSONDecodeError as e:
            return JsonResponse({'status': 'error', 'message': f'Error decoding JSON in upload_matrix view: {str(e)}'})
    else:
        return JsonResponse({'status': 'error', 'message': 'Invalid request method in upload_matrix view'})
    

    
