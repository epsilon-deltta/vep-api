import sys, types, string
import random
import re

bases = ['A','T','C','G']
baseH = dict([(x,i) for i,x in enumerate(bases)])
baseF = lambda x: baseH[x]
dna2fourthdecimal= lambda x: ''.join(map(str,map(baseF,list(x))))
dna2decimal = lambda x: int(dna2fourthdecimal(x),4)

wc_DNA = {'N':'N','.':'.','C':'G','G':'C','A':'T','T':'A','*':'*'}
wc_RNA = {'N':'N','.':'.','C':'G','G':'C','A':'U','U':'A','*':'*'}

def diff(templ,query):

    if len(templ)!=len(query):
        print(templ,query)
        raise Exception

    result = []
    nmatch = 0

    for i in range(len(query)):
        if templ[i]==query[i]:
            result.append(query[i])
            nmatch+=1
        else:
            result.append(query[i].lower())

    return ''.join(result),nmatch

def rev(seq):

#    if not isinstance(seq,types.StringType):
#        print "error: string type required!"
#        print type(seq)
#        sys.exit(1)

    seq_rev = ""

    for i in range(0,len(seq)):
        seq_rev = seq[i] + seq_rev

    return seq_rev


def compl(seq, seqType='DNA'):

    if not isinstance(seq,str):
        print ("error in mybio.compl: string type required!")
        print ("you entered:%s" % (seq,))
        print (type(seq))
        sys.exit(1)

    seq_ret = []

    for c in seq.upper():
        if seqType=='DNA':
            try:
                tmp = wc_DNA[c]
            except KeyError:
                tmp = c
            seq_ret.append(tmp)
        elif seqType=='RNA':
            try:
                tmp = wc_RNA[c]
            except KeyError:
                tmp = c
            seq_ret.append(tmp)
        else:
            raise Exception

    return ''.join(seq_ret)

def rc(seq,seqType='DNA'):

#    if not isinstance(seq,types.StringType):
#        print "error: string type required!"
#        print type(seq)
#        sys.exit(1)

    seq = rev(seq)
    seq = compl(seq,seqType)

    return seq


def motif_to_regex(motif):

    return motif.replace('U','T').replace('W','[AT]').replace('S','[CG]').replace('M','[AC]').replace('K','[GT]').replace('R','[AG]').replace('Y','[CT]').replace('B','[CGT]').replace('D','[AGT]').replace('H','[ACT]').replace('V','[ACG]').replace('N','[ACGT]') 


def count_regex(regex, seq):

    return sum(1 for _ in re.finditer('(?=(%s))' % regex, seq))

def random_dna_list(length,count=1,rna=False):

    seqL = []

    for i in range(count):
        if rna:
            seqL.append(''.join(random.choice('CGUA') for _ in xrange(length)))
        else:
            seqL.append(''.join(random.choice('CGTA') for _ in xrange(length)))

    return seqL

def decision(probability):
    return random.random() < probability

def spike_in_motif(seqL,motif,freq):

    resultL = []

    for s in seqL:
        
        while motif in s:
            s = s.replace(motif,''.join(random.choice('CGTA') for _ in xrange(len(motif))))

        if decision(freq):
            s_new = bytearray(s)
            posSta = random.randint(0,len(s)-len(motif))
            s_new[posSta:posSta+len(motif)] = motif
            s = str(s_new)

        resultL.append(s)

    return resultL


def switch_base(seq,idx,ref,alt):

    seq_new = bytearray(seq)

    if seq[idx] != ref:
        print (seq,idx,ref,alt)
        raise Exception

    seq_new[idx] = alt

    return str(seq_new)

def composition(seq,unit='percent'):

    import collections

    cntr = collections.Counter(seq)
    h = {}

    if unit=='percent':
        for k,v in cntr.items():
            h[k] = '%.3f' % (float(v)/sum(cntr.values()),)
    else:
        h = dict(cntr)

    return h
