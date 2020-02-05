from django.urls import path
import expression_profiles.views as views

app_name = 'expression_profiles'
urlpatterns = [
    path('expression_profiles/', views.expression_profiles, name='expression_profiles')
]
