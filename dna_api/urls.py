from django.urls import path,include

from . import views
# app_name = 'dna_api'
urlpatterns = [
    path('', views.index, name='index'),
    # path('<int:question_id>/', views.detail, name='detail'),
    # # ex: /polls/5/results/
    # path('<int:question_id>/results/', views.results, name='results'),
    # path('<int:question_id>/vote/', views.vote, name='vote'),
]