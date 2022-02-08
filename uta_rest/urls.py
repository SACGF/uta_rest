"""
The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.0/topics/http/urls/
"""
from django.urls import path

from uta_rest import views

urlpatterns = [
    path('', views.index, name='index'),
    path('transcript/<transcript_version>', views.transcript, name='transcript'),
    path('transcript/<uta_version>/<transcript_version>', views.uta_transcript, name='uta_transcript'),
]
