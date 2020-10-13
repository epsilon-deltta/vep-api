from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import get_object_or_404,render
from django.http import HttpResponse ,HttpRequest
import simplejson 
import dna.dnaMaker as maker 
import json
import dna.seq_processor as sp

# from rest_framework.views import APIView
# from rest_framework.response import Response
# from .serializers import UserSerializer
# from rest_framework import status


# Create your views here.

def shaodi(request):
    # latest_question_list = Question.objects.order_by('-pub_date')[:5]
    # context = {'latest_question_list': latest_question_list}
    return render(request, 'dna_api/index.html', context)
    # response = HttpResponse('success')
    # response["Access-Control-Allow-Methods"] = "GET, POST, OPTIONS"
    # response["Access-Control-Allow-Origin"] = "*"
    # response["Acess-Control-Max-Age"] = "1000"
    # response["Access-Control-Allow-Headers"] = "X-Requested-With, Content-Type"
#  seq/test
# @csrf_exempt 
def index(request):
    # latest_question_list = Question.objects.order_by('-pub_date')[:5]
    # context = {'latest_question_list': latest_question_list}
    context        = {}
    dna            = maker.get_random_dna()
    context['dna'] = dna 
    context.update(  maker.get_percentage(dna) )
    
    response        = HttpResponse(simplejson.dumps(context), content_type="application/json")

    return  response
    # response = HttpResponse('success')
    # response["Access-Control-Allow-Methods"] = "GET, POST, OPTIONS"
    # response["Access-Control-Allow-Origin"] = "*"
    # response["Acess-Control-Max-Age"] = "1000"
    # response["Access-Control-Allow-Headers"] = "X-Requested-With, Content-Type"
    # latest_question_list = Question.objects.order_by('-pub_date')[:5]
    # context = {'latest_question_list': latest_question_list}
@csrf_exempt 
#  seq/    
def seq(request):
    context  = {}
    if request.body is not None:
        data = json.loads(request.body)
        context  = sp.get_hexamer_track_info(data['seq'])
    else:
        print("It's empty")
    response = HttpResponse(simplejson.dumps(context), content_type="application/json")
    return response

@csrf_exempt 
def test_client_json(request):
    print(request.body)
    data = json.loads(request.body)
    print(data['test'])
    if request.is_ajax():
        print("==========")
        if request.method == 'POST':
            print('Raw Data: %s' % request.body)   
    return HttpResponse("OK")
# /seq/spliceai
@csrf_exempt 
def spliceAi(request):
    context ={}
    if request.body is not None:
        data = json.loads(request.body)
        context  = sp.get_splice_ai(data['seq'])
    else:
        print("It's empty")
    response = HttpResponse(simplejson.dumps(context), content_type="application/json")
    return response

# /seq/spliceai_opt
@csrf_exempt 
def spliceAi_option(request):
    context ={}
    if request.body is not None:
        data = json.loads(request.body)
        context  = sp.get_splice_ai_option(data['seq'],data['option'])
    else:
        print("It's empty")
    response = HttpResponse(simplejson.dumps(context), content_type="application/json")
    return response
    
# /seq/hex_mas
@csrf_exempt 
def hex_mas(request):
    context ={}
    if request.body is not None:
        data = json.loads(request.body)
        context  = sp.get_hexamer_maxent1(data['seq'])
    else:
        print("It's empty")
    response = HttpResponse(simplejson.dumps(context), content_type="application/json")
    return response
# /seq/test0
def test(request):
    return HttpResponse("Test successed")
    