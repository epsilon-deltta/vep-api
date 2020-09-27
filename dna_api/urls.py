from django.urls import path,include

from . import views
# app_name = 'dna_api'
urlpatterns = [
    path('', views.seq, name='seq'),
    path('test', views.index, name='index'),
    path('client', views.test_client_json, name='index'),
    path('splice',views.spliceAi,name="spliceAi"),
    # path('<int:question_id>/', views.detail, name='detail'),
    # # ex: /polls/5/results/
    # path('<int:question_id>/results/', views.results, name='results'),
    # path('<int:question_id>/vote/', views.vote, name='vote'),
    # path('', views.seq.as_view(), name='seq'),
]