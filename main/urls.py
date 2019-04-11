from django.urls import path
import main.views as main_views

urlpatterns = [
        path('', main_views.home),
        path('user_guide/', main_views.user_guide),
        path('contact/', main_views.contact)
]
