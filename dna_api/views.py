from django.shortcuts import get_object_or_404,render
from django.http import HttpResponse ,HttpRequest
import simplejson as json
import dna.dnaMaker as maker 
# from .models import Question
# Create your views here.

def shaodi(request):
    # latest_question_list = Question.objects.order_by('-pub_date')[:5]
    # context = {'latest_question_list': latest_question_list}
    return render(request, 'dna_api/index.html', context)

def index(request):
    # latest_question_list = Question.objects.order_by('-pub_date')[:5]
    # context = {'latest_question_list': latest_question_list}
    context        = {}
    dna            = maker.get_random_dna()
    context['dna'] = dna 
    context.update(  maker.get_percentage(dna) )
    

    return HttpResponse(json.dumps(context), content_type="application/json")