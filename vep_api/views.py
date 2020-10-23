from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import get_object_or_404,render
from django.http import HttpResponse ,HttpRequest
import simplejson 
import dna.dnaMaker as maker 
import json
from . import vep_request as vr
# import dna.seq_processor as sp

# from rest_framework.views import APIView
# from rest_framework.response import Response
# from .serializers import UserSerializer
# from rest_framework import status


# Create your views here.


    # response = HttpResponse('success')
    # response["Access-Control-Allow-Methods"] = "GET, POST, OPTIONS"
    # response["Access-Control-Allow-Origin"] = "*"
    # response["Acess-Control-Max-Age"] = "1000"
    # response["Access-Control-Allow-Headers"] = "X-Requested-With, Content-Type"
    # latest_question_list = Question.objects.order_by('-pub_date')[:5]
    # context = {'latest_question_list': latest_question_list}
@csrf_exempt 
#  seq/    
# def seq(request):
#     context  = {}
#     if request.body is not None:
#         data = json.loads(request.body)
#         context  = sp.get_hexamer_track_info(data['seq'])
#     else:
#         print("It's empty")
#     response = HttpResponse(simplejson.dumps(context), content_type="application/json")
#     return response

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
# /test
def test(request):
    return HttpResponse("Test successed")
#output:[{...},{...},..] {...} = Uploaded_variation Location Allele Gene Feature Feature_type Consequence cDNA_position CDS_position Protein_position Amino_acids Codons Existing_variation Extra
@csrf_exempt 
def everything(request):
    context ={}
    if request.body is not None:
        data = json.loads(request.body)
        context  = vr.get_vep_data(data['seq'])
    else:
        print("It's empty")
    response = HttpResponse(simplejson.dumps(context), content_type="application/json")
    return response
    
# @csrf_exempt 
# def sift_poly(request):
#     context ={}
#     if request.body is not None:
#         data = json.loads(request.body)
#         context  = vr.get_vep_data(data['seq'],option='-sift b --polyphen p ')
#     else:
#         print("It's empty")
#     response = HttpResponse(simplejson.dumps(context), content_type="application/json")
#     return response