from django.shortcuts import render, render_to_response
from main import models as main_models

# Create your views here.

def home(request):
    news = main_models.News.objects.all().order_by('-release_date')
    return render_to_response('main/home.html',locals())

def user_guide(request):
    return render_to_response('main/user_guide.html',locals())

def contact(request):
    return render_to_response('main/contact.html',locals())
