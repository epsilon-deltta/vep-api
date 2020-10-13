from django.urls import path,include

from . import views
# app_name = 'dna_api'
#/seq/~
urlpatterns = [
    path(''            , views.seq             , name='seq'),
    path('test'        , views.index           , name='index'),
    path('client'      , views.test_client_json, name='index'),
    path('spliceai'    , views.spliceAi        , name="spliceAi"),
    path('spliceai_opt', views.spliceAi_option , name="spliceAi_option"),
    path('hex_mas'     , views.hex_mas         , name="hex_mas"),
    path('test0',views.test,name="test0"),
    # path('<int:question_id>/', views.detail, name='detail'),
    # # ex: /polls/5/results/
    # path('<int:question_id>/results/', views.results, name='results'),
    # path('<int:question_id>/vote/', views.vote, name='vote'),
    # path('', views.seq.as_view(), name='seq'),
]