import sys, os, copy, re, gzip, collections
import pandas as pd
import numpy as np
import jkbasic, jkbio
# import socket
import random as rd

###############################################################
###############################################################
assemblyH = {'hg18':'/data1/Sequence/ucsc_hg18/hg18_nh.fa', 'hg19':'/data1/Sequence/ucsc_hg19/hg19_nh.fa', 'hg38':'/Z/DB/Homo_sapiens.GRCh38.dna.primary_assembly.fa'}

sep = os.path.sep
homedir=os.path.dirname(__file__) +sep + 'data'
refFlat_path = '%s/D/refFlat.txt' %(homedir) #added 
###############################################################
###############################################################
# ===============locus================================

# ===========

# ================splice AI=================== #
def spliceAI(locusStr): # locusStr: 1-base, 1-base

    result = []
    
    for r in spliceAI_raw(locusStr):
        tokL = r[:-1].split('\t')

        chrN = tokL[0]
        pos = int(tokL[1])

        ref = tokL[3] 
        alt = tokL[4]

        tL = tokL[-1].split('|')

        geneName = tL[1]

        score  = dict(zip(('AG','AL','DG','DL'),map(float,tL[2:6]) )  )
        relPos = dict(zip(('AG','AL','DG','DL'),map(int,tL[6:])    )  )

        result.append({'chrN':chrN,'pos':pos,'ref':ref,'alt':alt,'geneName':geneName,'score':score,'relPos':relPos})
    return result
def spliceAI_raw(locusStr): # locusStr: 1-base, 1-base

    if locusStr[:3] == 'chr':
        locusStr = locusStr[3:]
    #posSta,posEnd = map(int,locusStr.split(':')[-1].split('-'))

    return tbi_bed_query('%s/BigFiles/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz' % homedir, locusStr) + tbi_bed_query('%s/BigFiles/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz' % homedir, locusStr)


def tbi_bed_query(bedFilePath,locusStr):
    if os.path.isfile(bedFilePath) == False:
        print("this is not file | There is not file")
        raise FileNotFoundError
    # print('%s/tools/tabix-0.2.6/tabix %s %s' % (homedir,bedFilePath,locusStr))    
    f = os.popen('%s/tools/tabix-0.2.6/tabix %s %s' % (homedir,bedFilePath,locusStr),'r')
    result = f.readlines()
    # print("result" ,result)
    return result

def primateAI_raw(locusStr): # locusStr: 1-base, 1-base

    if locusStr[:3] != 'chr':
        locusStr = 'chr'+locusStr

    chrom = locusStr.split(':')[0]
    posSta, posEnd = map(int,locusStr.split(':')[-1].split('-'))

    locusStr_new = '%s:%s-%s' % (chrom,posSta-1,posSta+1)

    return tbi_bed_query('%s/BigFiles/primateAI/primateAI_hg38.bed.gz' % homedir, locusStr_new)

###############################################################
###############################################################
# =================hexamer_track==================== #
def hexamer_track(seq,hexH4=None,figPath='default'):

    if not hexH4:
        hexH4 = load_hexamer4()
 
    from matplotlib.font_manager import FontProperties
    import pylab

    #pylab.rc('mathtext', default='regular')

    font = FontProperties()
    font.set_family('monospace')
    font.set_size(17.*84./50*7/15)       

    fig = pylab.figure(figsize=(7*len(seq)/50.,7))

    regionType = ['5p_exon','5p_intron','3p_exon','3p_intron']

    scoreH = {}

    for r in range(4):

        ax = fig.add_subplot('41%s' % (r+1,))
        ax.axis([0,len(seq),0,6])
        ax.set_ylabel(regionType[r])
        
        if r==0:
            ax.set_title('Hexamer score')
        
        scoreH[regionType[r]] = {}

        for i in range(len(seq)-5):
            s = seq[i:i+6]
            c = calColor(hexH4[regionType[r]][s])
            scoreH[regionType[r]][s] = hexH4[regionType[r]][s]
            ax.text(i,i%6+0.2,s,color=c,fontproperties=font)
        
    if figPath == 'default':
        num = rd.randrange(0,500)   # it's for rainy day ,if many people access file,..
        sp=os.path.sep
        savedPath ='%shexamer_track'% (os.path.dirname(__file__)+sp+'output'+sp)+str(num)+'.png' 
        print("saved path :",savedPath)
        
        fig.savefig(savedPath)
        # fig.savefig('%s/hexamer_track.pdf' % homedir) //original
        scoreH['path']=savedPath 
    elif figPath:
        fig.savefig(figPath)
        scoreH['path']=figPath
    return scoreH

def load_hexamer4():

    result = {}

    for h in ['5p_exon','3p_exon','5p_intron','3p_intron']:
        result[h] = load_hexamer(h)
    
    return result

def load_hexamer(segType):

    if segType == '5p_exon':
        s = ('A5','R1')
    elif segType == '5p_intron':
        s = ('A5','R2')
    elif segType == '3p_intron':
        s = ('A3','R1')
    elif segType == '3p_exon':
        s = ('A3','R2')
    else:
        print('segType unknown: %s\n' % segType)
        raise Exception

    import dnatools
    sep = os.path.sep #ypil added
    return pd.Series(index=dnatools.make_mer_list(6),data=np.load('%sPackages/cell-2015/results/N4_Motif_Effect_Sizes/%sSS/mean_effects_sdpos_%s.npy' % (os.path.dirname(__file__)+sep+'data'+sep,s[0],s[1])))
# ==============================================================
def hexamer4_byCoord_general(chrNum,chrSta,chrEnd,strand,ref,alt,hexH4,assembly): # pos, 0-base, inclusive

    ref = '' if ref=='-' else ref
    alt = '' if alt=='-' else alt

    l = locus('chr%s:%s-%s%s' % (chrNum,chrSta-6,chrEnd+5,strand))

    s_before = l.twoBitFrag(assembly)

    if s_before[5:5+len(ref)] != ref:
        print(s_before, ref, alt)
        raise Exception

    s_after = s_before[:5] + alt + s_before[-5:]

    if 'N' in s_before or 'N' in s_after:
        print (s_before,s_after)
        raise Exception

    result = {}

    for h,hexH in hexH4.items():
        
        scores_after = hexamer(s_after,hexH)
        scores_before = hexamer(s_before,hexH)

        result[h] = sum(scores_after) - sum(scores_before)

    #return scores_before, scores_after, scores_after - scores_before, sum(scores_after - scores_before)
    return result

def hexamer(seq,hexH): 

    resultL = []

    for i in range(len(seq)-5):
        
        s = seq[i:i+6]
        resultL.append(hexH[s])

    return pd.Series(resultL)
###############################################################
###############################################################
# input: 12,51391349,'T','GG'
# output_exam 
# { 'count': 2, 
#      0   : { 'type'     : 'snv', 
#              'strand'   : '-'    ,         
#              'transName': 'GALNT6',
#              'hexamer'  : { '5p_exon'  : 0.26357736584766656, 
#                             '3p_exon'  : -1.519055738781958, 
#                             '5p_intron': -0.016999846994960266, 
#                             '3p_intron': -1.7098262852115953 }, 
#              'mes'      : { 'donor'   : {'delta': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
#                                          'final': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, 
#                             'acceptor': { 'delta': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
#                                           'final': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}}             },
#      1   : { ...
#              'mes'      : { 'donor'   : {'ref': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
#                                          'alt': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, 
#                             'acceptor': { 'ref': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
#                                           'alt': [0.0, 0.0, 0.0, 0. ....
# }


def variant_bi2(chrNum,chrPos,ref,alt,assembly='hg38',hexH4=False):
    # ref>alt conversion: + strand oriented
    
    # ref alt: unifying nomenclature
    uni_nomen = lambda x: x.replace('-','').replace('.','').replace('*','')
    ref = uni_nomen(ref)
    alt = uni_nomen(alt)
    
    if not hexH4:
        hexH4 = load_hexamer4()
    
    regions = locus("%s:%s-%s%s"%('chr'+str(chrNum), chrPos-1, chrPos, '+')).regionType() # UCSC genome 0-based
    #(t['transName'],t['transID'],flag, sense, marg)
    r_ls = []
    for region in regions:
        s = '+' if region[3] =='sense' else '-'
        if region[2]!='lnc' and (region[0],s) not in r_ls:
            r_ls.append((region[0],s))
            
    r_ls = list(set(r_ls))
    
    result ={}
    result['count'] = len(r_ls)
    
    if not len(r_ls): # only 'lnc' regions <= how to deal with?
        print('Given region has long non coding (lnc) annotations only')
    for i,r in enumerate(r_ls):
        result[i] = {}

        transName, strand = r

        result[i]['transName'] = transName
        result[i]['strand']    = strand

        if strand == '-':
            ref_, alt_ = jkbio.rc(ref), jkbio.rc(alt)
        else:
            ref_, alt_ = ref, alt
        
        hex_result = hexamer4_byCoord_general(chrNum,chrPos,chrPos,strand,ref_,alt_,hexH4,assembly)
        result[i]['hexamer'] = hex_result

        # result['hexamer'] = hex_result
        if len(ref)==len(alt):

            result[i]['type']    = 'snv'
            
            mes_result = mes_byCoord(chrNum,chrPos,strand,ref_,alt_,assembly,verbose=True)
            # 'MES donor'
            donor_delta = list(map(float,list(mes_result[0])))
            donor_final = list(map(float,mes_result[2]))
            # 'MES acceptor'
            acceptor_delta = list(map(float,list(mes_result[1])))
            acceptor_final = list(map(float,mes_result[3]))

            maxents_dict = {}
            maxents_dict['donor'   ] = {'delta' :donor_delta,'final' : donor_final }   
            maxents_dict['acceptor'] = {'delta' :acceptor_delta,'final' : acceptor_final }
            result[i]['mes'] = maxents_dict 
        else:

            result[i]['type'] = 'indel'

            mes_result = mes_byCoord_general(chrNum,chrPos,chrPos+len(ref_)-1,strand,ref_,alt_,assembly)
            # 'MES donor' 
            donor_ref = list(map(float,mes_result[0][0]))
            donor_alt = list(map(float,mes_result[0][1]))
            # 'MES acceptor'
            acceptor_ref = list(map(float,mes_result[1][0]))
            acceptor_alt = list(map(float,mes_result[1][1]))

            
            maxents_dict = {}
            maxents_dict['donor'   ] = {'ref' :donor_ref,'alt' : donor_alt }   
            maxents_dict['acceptor'] = {'ref' :acceptor_ref,'alt' : acceptor_alt }
            result[i]['mes'] = maxents_dict
    return result 
# ========================================================
# ==================maxentscan============================
def mes_byCoord(chrNum,pos,strand,ref,alt,assembly,verbose=False): # assuming single nt to single nt change, base-1 inclusive

    l = locus('chr%s:%s-%s%s' % (chrNum,pos-23,pos+22,strand))

    s_before = l.twoBitFrag(assembly)

    if s_before[22:22+len(ref)] != ref:
        print(s_before, ref, alt)
        raise Exception

    s_after = s_before[:22] + alt + s_before[-22:]

    mes5_tup = (mes5_scan(s_before[14:-14]),mes5_scan(s_after[14:-14]),)
    mes3_tup = (mes3_scan(s_before),mes3_scan(s_after),)

    mes5_dt = np.array(mes5_tup[1])-np.array(mes5_tup[0])
    mes3_dt = np.array(mes3_tup[1])-np.array(mes3_tup[0])

#    if verbose==True:
#        print(mes5_tup,mes3_tup)

    if verbose:
        return mes5_dt,mes3_dt,mes5_tup[1],mes3_tup[1]
    else:
        return mes5_dt,mes3_dt

def mes_byCoord_general(chrNum,chrSta,chrEnd,strand,ref,alt,assembly): # base-1 inclusive

    ref = '' if ref=='-' else ref
    alt = '' if alt=='-' else alt

    l = locus('chr%s:%s-%s%s' % (chrNum,chrSta-23,chrEnd+22,strand))

    s_before = l.twoBitFrag(assembly)

    if s_before[22:22+len(ref)] != ref:
        print (s_before, chrNum, chrSta, chrEnd, ref, alt)
        raise Exception

    s_after = s_before[:22] + alt + s_before[-22:]

    mes5_tup = (mes5_scan(s_before[14:-14]),mes5_scan(s_after[14:-14]))
    mes3_tup = (mes3_scan(s_before),mes3_scan(s_after))

    return mes5_tup,mes3_tup
# def mes_byCoord_general(chrNum,chrSta,chrEnd,strand,ref,alt,assembly): # base-0 exclusive

#     ref = '' if ref=='-' else ref
#     alt = '' if alt=='-' else alt

#     l = locus('chr%s:%s-%s%s' % (chrNum,chrSta-22,chrEnd+22,strand))

#     s_before = l.twoBitFrag(assembly)

#     if s_before[22:22+len(ref)] != ref:
#         print(s_before, chrNum, chrSta, chrEnd, ref, alt)
#         raise Exception

#     s_after = s_before[:22] + alt + s_before[-22:]

#     mes5_tup = (mes5_scan(s_before[14:-14]),mes5_scan(s_after[14:-14]),)
#     mes3_tup = (mes3_scan(s_before),mes3_scan(s_after))

#     return mes5_tup,mes3_tup

def mes_byCoord_general_b1(chrNum,chrSta,chrEnd,strand,ref,alt,assembly): # base-0 exclusive
    chrSta = chrSta -1
    chrEnd = chrEnd -2

    ref = '' if ref=='-' else ref
    alt = '' if alt=='-' else alt

    l = locus('chr%s:%s-%s%s' % (chrNum,chrSta-22,chrEnd+22,strand))

    s_before = l.twoBitFrag(assembly)

    if s_before[22:22+len(ref)] != ref:
        print(s_before, chrNum, chrSta, chrEnd, ref, alt)
        raise Exception

    s_after = s_before[:22] + alt + s_before[-22:]

    mes5_tup = (mes5_scan(s_before[14:-14]),mes5_scan(s_after[14:-14]),)
    mes3_tup = (mes3_scan(s_before),mes3_scan(s_after))

    return mes5_tup,mes3_tup    

def mes5_scan(seq):

    seq = seq.upper()

    resultL = [None]*3

    for i in range(len(seq)-8):
        if seq[i+3:i+5] == 'GT':
            resultL.append(mes5(seq[i:i+9]))
        else:
            resultL.append(None)

    resultL += [None]*5

    resultL = [(0 if x==None or x<0 else x) for x in resultL]

    return resultL

def mes3_scan(seq):

    seq = seq.upper()

    resultL = [None]*19

    for i in range(len(seq)-22):
        if seq[i+18:i+20] == 'AG':
            resultL.append(mes3(seq[i:i+23]))
        else:
            resultL.append(None)

    resultL += [None]*3

    resultL = [(0 if x==None or x<0 else x) for x in resultL]

    return np.array(resultL)
def mes5(seq): # seq example: CAGgtaagt

    if len(seq) != 9:
        raise Exception
    print(os.popen('cd %s/tools/maxentscan; perl score5_mod.pl %s' % (homedir,seq)).readline())    
    return float(os.popen('cd %s/tools/maxentscan; perl score5_mod.pl %s' % (homedir,seq)).readline())

def mes3(seq): # seq example: ttcca aacga acttt tgtag GGA (23)

    if len(seq) != 23:
        raise Exception

    return float(os.popen('cd %s/tools/maxentscan; perl score3_mod.pl %s' % (homedir,seq)).readline())




def mes_byCoord_subs(chrNum,chrSta,chrEnd,strand,ref,alt,assembly,verbose=False): # assuming substitution, base-1 inclusive

    l = locus('chr%s:%s-%s%s' % (chrNum,chrSta-23,chrEnd+22,strand))

    s_before = l.twoBitFrag(assembly)

    if s_before[22:22+len(ref)] != ref:
        print(s_before, ref, alt)
        raise Exception

    s_after = s_before[:22] + alt + s_before[-22:]

    mes5_tup = (mes5_scan(s_before[14:-14]),mes5_scan(s_after[14:-14]),)
    mes3_tup = (mes3_scan(s_before),mes3_scan(s_after),)

    mes5_dt = np.array(mes5_tup[1])-np.array(mes5_tup[0])
    mes3_dt = np.array(mes3_tup[1])-np.array(mes3_tup[0])

#    if verbose==True:
#        print(mes5_tup,mes3_tup)

    if verbose:
        return mes5_dt,mes3_dt,mes5_tup,mes3_tup
    else:
        return mes5_dt,mes3_dt


def mes5_byCoord(chrNum,pos,strand,assembly): # pos is 1-base sta or end of intron

    if strand == '+':
        l = locus('chr%s:%s-%s%s' % (chrNum,pos-4,pos+5,'+'))
    else:
        l = locus('chr%s:%s-%s%s' % (chrNum,pos-6,pos+3,'-'))

    s = l.twoBitFrag(assembly)
    return s,mes5(s)

def mes3_byCoord(chrNum,pos,strand,assembly): # pos is 1-base sta or end of intron

    if strand == '+':
        l = locus('chr%s:%s-%s%s' % (chrNum,pos-20,pos+3,'+'))
    else:
        l = locus('chr%s:%s-%s%s' % (chrNum,pos-4,pos+19,'-'))

    s = l.twoBitFrag(assembly)
    return s,mes3(s)
# ======================
# ======================
def loadFasta(fastaPath='%s/D/Sequences/gencode_24/gencode.v24lift37.pc_transcripts.fa.gz' % homedir, blacklist=['NR_106988']):

    h = collections.defaultdict(list)    

    if fastaPath.endswith('.gz'):
        f = gzip.open(fastaPath)
    else:
        f = open(fastaPath)

    seqID = None
    seqL = []

    while 1:

        line = f.readline()

        if not line or line.startswith('>'):

            if seqID and seqID not in blacklist:
                h[seqID].append(''.join(seqL))
                seqL = []

            if line:
                seqID = line.split('|')[0].split(' ')[0][1:]
            else:
                break

        else:

            seqL.append(line.rstrip())

    return h

def loadBlatOutput(blatOutputPath,by='transID',blacklist=['NR_106988']):

    h = collections.defaultdict(list)
    if blatOutputPath.endswith('.gz'):
        f = gzip.open(blatOutputPath)
    else:
        f = open(blatOutputPath)

    for line in f:
        
        if line[0] == '#':
            continue

        r = processBlatLine(line)

        if r['transID'] in blacklist:
            continue

        h[r[by]].append(r)

    from operator import attrgetter, itemgetter

    for k,vL in list(h.items()):
        h[k] = sorted(vL,key=itemgetter('txnSta','txnEnd'))
    return h
def loadAppris_refseq(path = '%s/D/Sequences/APPRIS_RefSeq107_20170423_src20170406' % homedir):

    h = collections.defaultdict(list)

    for line in open(path):

        tokL = line[:-1].split('\t')

        if tokL[2][:3] == 'NM_' and tokL[4].startswith('PRINCIPAL'):
            h[tokL[0]].append((tokL[2],tokL[4]))

    for geneN,transL in h.items():

        transL.sort(lambda x,y: cmp(x[1],y[1]))
        transL.sort(lambda x,y: cmp(int(x[0][3:].split('.')[0]),int(y[0][3:].split('.')[0])))
        h[geneN] = transL[0][0]

    return h

# gencode basic: D/Sequences/gencode_24/wgEncodeGencodeBasicV24lift37.txt.gz
# gencode comp: D/Sequences/gencode_24/wgEncodeGencodeCompV24lift37.txt.gz

# def loadBlatOutput(blatOutputPath,by='transID',blacklist=['NR_106988']):

#     h = collections.defaultdict(list)
#     if blatOutputPath.endswith('.gz'):
#         f = gzip.open(blatOutputPath)
#     else:
#         f = open(blatOutputPath)

#     for line in f:

#         if line[0] == '#':
#             continue
    
#         r = processBlatLine(line)

#         if r['transID'] in blacklist:
#             continue

#         h[r[by]].append(r)

#     from operator import itemgetter, attrgetter

#     for k,vL in list(h.items()):
#         h[k] = sorted(vL,key=itemgetter('txnSta','txnEnd'))
#     return h

def loadBlatOutputByGene(blatOutputPath=refFlat_path):

    return loadBlatOutput(blatOutputPath, 'geneName')

def loadBlatOutputByID(blatOutputPath=refFlat_path):

    return loadBlatOutput(blatOutputPath,'transID')
# conver genome to trans
def loadBlatOutputByChr(blatOutputPath=refFlat_path):

    return loadBlatOutput(blatOutputPath, 'chrom')
    
def processBlatLine(line):

    tokL = line.rstrip().split('\t')

    h = {}

    h['transName'] = tokL[0]
    h['transID'] = tokL[1]
    h['chrom'] = tokL[2]
    h['chrNum'] = tokL[2][3:]
    h['strand'] = tokL[3]
    h['txnSta'] = int(tokL[4])
    h['txnEnd'] = int(tokL[5])
    h['txnLen'] = h['txnEnd'] - h['txnSta']
    h['cdsSta'] = int(tokL[6])
    h['cdsEnd'] = int(tokL[7])
    h['exnList'] = list(map(lambda x,y: (int(x),int(y)), tokL[9].split(',')[:-1], tokL[10].split(',')[:-1]))
    h['exnLenList'] = [e-s for (s,e) in h['exnList']]
    h['exnLen'] = sum(h['exnLenList'])

    if len(tokL) > 12:
        h['geneName'] = tokL[12]

    h['cdsList'] = []
    frontL, backL = [],[]

    if h['cdsSta'] != h['cdsEnd']:

        for (s,e) in h['exnList']:

            frontL.append((min(s,h['cdsSta']),min(e,h['cdsSta'])))
            h['cdsList'].append((max(s,h['cdsSta']),min(e,h['cdsEnd'])))
            backL.append((max(s, h['cdsEnd']), max(e, h['cdsEnd'])))

        frontL = [x for x in frontL if x[0] < x[1]]
        h['cdsList'] = [x for x in h['cdsList'] if x[0] < x[1]]
        backL = [x for x in backL if x[0] < x[1]]

    if h['strand'] == '+':
            h['utr5'] = frontL
            h['utr3'] = backL
    elif h['strand'] == '-':
            h['utr5'] = backL
            h['utr3'] = frontL
    else:
            raise Exception

    h['utr5Len'] = sum([e-s for (s,e) in h['utr5']])
    h['utr3Len'] = sum([e-s for (s,e) in h['utr3']])
    h['cdsLen'] = sum([e-s for (s,e) in h['cdsList']])

    h['intron'] = []

    for i in range(len(h['exnList'])-1):
        h['intron'].append((h['exnList'][i][1],h['exnList'][i+1][0]))

    return h
# old version
# def processBlatLine(line):

#     tokL = line.rstrip().split('\t')

#     h = {}

#     h['transName'] = tokL[0]
#     h['transID'] = tokL[1]
#     h['chrom'] = tokL[2]
#     h['chrNum'] = tokL[2][3:]
#     h['strand'] = tokL[3]
#     h['txnSta'] = int(tokL[4])
#     h['txnEnd'] = int(tokL[5])
#     h['txnLen'] = h['txnEnd'] - h['txnSta']
#     h['cdsSta'] = int(tokL[6])
#     h['cdsEnd'] = int(tokL[7])
#     h['exnList'] = list(map(lambda x,y: (int(x),int(y)), tokL[9].split(',')[:-1], tokL[10].split(',')[:-1]))
#     h['exnLenList'] = [e-s for (s,e) in h['exnList']]
#     h['exnLen'] = sum(h['exnLenList'])

#     if len(tokL) > 12:
#         h['geneName'] = tokL[12]

#     h['cdsList'] = []
#     frontL, backL = [],[]

#     if h['cdsSta'] != h['cdsEnd']:

#         for (s,e) in h['exnList']:

#             frontL.append((min(s,h['cdsSta']),min(e,h['cdsSta'])))
#             h['cdsList'].append((max(s,h['cdsSta']),min(e,h['cdsEnd'])))
#             backL.append((max(s, h['cdsEnd']), max(e, h['cdsEnd'])))

#         frontL = [x for x in frontL if x[0] < x[1]]
#         h['cdsList'] = [x for x in h['cdsList'] if x[0] < x[1]]
#         backL = [x for x in backL if x[0] < x[1]]
#         # frontL = filter(lambda x: x[0] < x[1],frontL)
#         # h['cdsList'] = filter(lambda x: x[0] < x[1], h['cdsList'])
#         # backL = filter(lambda x: x[0] < x[1], backL)

#     if h['strand'] == '+':
#             h['utr5'] = frontL
#             h['utr3'] = backL
#     elif h['strand'] == '-':
#             h['utr5'] = backL
#             h['utr3'] = frontL
#     else:
#             raise Exception

#     h['utr5Len'] = sum([e-s for (s,e) in h['utr5']])
#     h['utr3Len'] = sum([e-s for (s,e) in h['utr3']])
#     h['cdsLen'] = sum([e-s for (s,e) in h['cdsList']])

#     h['intron'] = []

#     for i in range(len(h['exnList'])-1):
#         h['intron'].append((h['exnList'][i][1],h['exnList'][i+1][0]))

#     return h


def margin(innerRange,outerRange): # 0-base

    i_s,i_e = innerRange
    o_s,o_e = outerRange

    return i_s-o_s, o_e-i_e+1



def getRegionTypeUsingTransH(transH,gPos): # transH: individual transcript hash; gPos: 1-base

    if len(transH['cdsList']) == 0:
        return 'lnc'

    if sum(map(lambda e: int(e[0]<gPos<=e[1]),transH['utr3'])) > 0:
        return 'utr3'
    elif sum(map(lambda e: int(e[0]<gPos<=e[1]),transH['utr5'])) > 0:
        return 'utr5'
    elif sum(map(lambda e: int(e[0]<gPos<=e[1]),transH['cdsList'])) > 0:
        return 'cds'
    else:
        print (transH,gPos)
        raise Exception


# nm_~ ,t_position
def convertTrans2Genome(transID,transPos,blatH=loadBlatOutputByID(),transLen=-1): # pos base-1

    result = []

    for trans in blatH[transID]:

        if '_' in  trans['chrom']:
            continue

        if transLen == -1:
            transLen = trans['exnLen']

        if trans['strand'] == '-':
            transPos = transLen-transPos+1

        tally = 1

        for i,n in enumerate(trans['exnLenList']):

            if tally <= transPos < tally+n:
                result.append((trans['chrom'], trans['exnList'][i][0]+(transPos-tally)+1, trans['strand'])) # pos base-1
                break

            tally += n

    return result
# 11:1234
# transcript - position

def convertGenome2Trans(chrNum,gPos,blatH_byChr=loadBlatOutputByChr() ): # genomic pos base1, transcript pos base0

    if 'chr'+str(chrNum) not in blatH_byChr:
        return []

    result = []

    for trans in blatH_byChr['chr'+str(chrNum)]:

        if trans['txnEnd'] < gPos:
            continue
        elif gPos <= trans['txnSta']:
            break

        for i,n in enumerate(trans['exnList']):

            if n[1] < gPos:
                continue
            elif gPos <= n[0]:
                break

            annot = getRegionTypeUsingTransH(trans, gPos)

            if trans['strand'] == '+':
                result.append((trans, sum([min(e[1],gPos)-min(e[0],gPos) for e in trans['exnList']])-1, annot))
            else:
                result.append((trans, sum([max(e[1],gPos)-max(e[0],gPos) for e in trans['exnList']]), annot))

            break

    return result

def loadRefFlatByChr(refFlatFileName='/%s/D/Sequences/hg19/refFlat_hg19.txt' % (homedir,)):

    h = {}

    for line in open(refFlatFileName):
    
        r = processBlatLine(line)
        jkbasic.addHash(h, r['chrom'], r)

    return h

def loadRefFlatByGeneName(refFlatFileName='/%s/D/Sequences/hg19/refFlat_hg19.txt' % (homedir,)):

    h = {}

    for line in open(refFlatFileName):
    
        r = processBlatLine(line)
        jkbasic.addHash(h, r['transName'], r)

    return h

# def processRefFlatLine(line):
#
#     tokL = line.rstrip().split('\t')
#
#     h = {}
#
#     h['geneName'] = tokL[0]
#     h['refSeqId'] = tokL[1]
#     h['chrom'] = tokL[2]
#     h['chrNum'] = tokL[2][3:]
#     h['strand'] = tokL[3]
#     h['txnSta'] = int(tokL[4])
#     h['txnEnd'] = int(tokL[5])
#     h['txnLen'] = h['txnEnd'] - h['txnSta']
#     h['cdsSta'] = int(tokL[6])
#     h['cdsEnd'] = int(tokL[7])
#     h['exnList'] = map(lambda x,y: (int(x),int(y)), tokL[9].split(',')[:-1], tokL[10].split(',')[:-1])
#     h['exnLenList'] = [e-s for (s,e) in h['exnList']]
#     h['exnLen'] = sum(h['exnLenList'])
#
#     h['cdsList'] = []
#     frontL, backL = [],[]
#
#     for (s,e) in h['exnList']:
#
#             if e<=h['cdsSta']:
#                     frontL.append((s,e))
#             elif s<=h['cdsSta'] and h['cdsSta']<=e:
#                     frontL.append((s,h['cdsSta']))
#                     h['cdsList'].append((h['cdsSta'],e))
#             elif h['cdsSta']<=s and e<=h['cdsEnd']:
#                     h['cdsList'].append((s,e))
#             elif s<=h['cdsEnd'] and h['cdsEnd']<=e:
#                     h['cdsList'].append((s,h['cdsEnd']))
#                     backL.append((h['cdsEnd'],e))
#             elif h['cdsEnd']<=s:
#                     backL.append((s,e))
#             else:
#                     raise Exception
#
#     if h['strand'] == '+':
#             h['utr5'] = frontL
#             h['utr3'] = backL
#     elif h['strand'] == '-':
#             h['utr5'] = backL
#             h['utr3'] = frontL
#     else:
#             raise Exception
#
#     return h
#
#
# def processKgLine(line):
#
#     tokL = line.rstrip().split('\t')
#
#     h = {}
#
#     h['geneId'] = tokL[0]
#     h['chrom'] = tokL[1]
#     h['chrNum'] = tokL[1][3:]
#     h['strand'] = tokL[2]
#     h['txnSta'] = int(tokL[3])
#     h['txnEnd'] = int(tokL[4])
#     h['txnLen'] = h['txnEnd'] - h['txnSta']
#     h['cdsSta'] = int(tokL[5])
#     h['cdsEnd'] = int(tokL[6])
#     h['exnList'] = map(lambda x,y: (int(x),int(y)), tokL[8].split(',')[:-1], tokL[9].split(',')[:-1])
#     h['exnLenList'] = [e-s for (s,e) in h['exnList']]
#     h['exnLen'] = sum(h['exnLenList'])
#
#     h['cdsList'] = []
#     frontL, backL = [],[]
#
#     for (s,e) in h['exnList']:
#
#         if e<=h['cdsSta']:
#             frontL.append((s,e))
#         elif s<=h['cdsSta'] and h['cdsSta']<=e:
#             frontL.append((s,h['cdsSta']))
#             h['cdsList'].append((h['cdsSta'],e))
#         elif h['cdsSta']<=s and e<=h['cdsEnd']:
#             h['cdsList'].append((s,e))
#         elif s<=h['cdsEnd'] and h['cdsEnd']<=e:
#             h['cdsList'].append((s,h['cdsEnd']))
#             backL.append((h['cdsEnd'],e))
#         elif h['cdsEnd']<=s:
#             backL.append((s,e))
#         else:
#             raise Exception
#
#     if h['strand'] == '+':
#         h['utr5pLen'] = sum([e-s for (s,e) in frontL])
#     elif h['strand'] == '-':
#         h['utr5pLen'] = sum([e-s for (s,e) in backL])
#     else:
#         raise Exception
#
#     h['cdsLen'] = sum([e-s for (s,e) in h['cdsList']])
#
#     if h['strand'] == '+':
#         exnLenListH = h['exnLenList']
#     else:
#         exnLenListH = h['exnLenList'][::-1]
#
#     transOffset = h['utr5pLen'] * -1
#     h['frame'] = []
#
#     # index of h['frame'] is [exon number]-1
#
#     for i in range(len(exnLenListH)):
#
#         if 0 <= transOffset < h['cdsLen']:
#             frame5p = transOffset % 3
#         else:
#             frame5p = None
#
#         if 0 <= transOffset+exnLenListH[i]-1 < h['cdsLen']:
#             frame3p = (transOffset+exnLenListH[i]-1) % 3
#         else:
#             frame3p = None
#
#         h['frame'].append((frame5p, frame3p))
#
#         transOffset += exnLenListH[i]
#
#     return h

#def parse_vcf_info(info):
#
#        itemL = info.split(';')
#        datH = {}
#
#        for item in itemL:
#                tag = item.split('=')[0]
#                val = '='.join(item.split('=')[1:])
#                datH[tag] = val
#
#        return(datH)


def parse_vcf_info(info):

    itemL = info.split(';')
    datH = {}

    for item in itemL:

        tag = item.split('=')[0]
        val = '='.join(item.split('=')[1:])

        if tag == 'GENE':

            rm = re.match('([^_]*)_?(ENS.[0-9]{11})',val)

            if rm:
                geneName = rm.group(1)
                ens_transID = rm.group(2)
                datH[tag] = (geneName,ens_transID)
            else:
                datH[tag] = (val,'')

        else:
            datH[tag] = val

    return(datH)

def getGenePos(refFlatFile='/data1/Sequence/ucsc_hg19/annot/refFlat.txt', geneList=[]):
    inFile = open(refFlatFile, 'r')
    posH = {}
    for line in inFile:
        colL = line[:-1].split('\t')
        gene_sym = colL[0]
        chrom = colL[2]
        pos = int(colL[4])
        if geneList==[] or gene_sym in geneList:
            if gene_sym not in posH:
                posH[gene_sym] = {'chrom':chrom, 'pos':pos}
            elif posH[gene_sym]['pos'] > pos:
                posH[gene_sym]['pos'] = pos
    return posH

# def loadLincByChr(dataFileName='/Z/Sequence/ucsc_hg19/annot/lincRNAsTranscripts.txt',h={}):
#
#     for line in open(dataFileName):
#
#         r = processLincLine(line)
#
#         jkbasic.addHash(h, r['chrom'], r)
#
#     return h
#
#
# def processLincLine(line):
#
#     tokL = line.rstrip().split('\t')[1:]
#
#     h = {}
#
#     h['geneId'] = tokL[0]
#     h['chrom'] = tokL[1]
#     h['chrNum'] = tokL[1][3:]
#     h['strand'] = tokL[2]
#     h['txnSta'] = int(tokL[3])
#     h['txnEnd'] = int(tokL[4])
#     h['txnLen'] = h['txnEnd'] - h['txnSta']
#     h['cdsSta'] = int(tokL[5])
#     h['cdsEnd'] = int(tokL[6])
#     h['exnList'] = map(lambda x,y: (int(x),int(y)), tokL[8].split(',')[:-1], tokL[9].split(',')[:-1])
#     h['exnLenList'] = [e-s for (s,e) in h['exnList']]
#     h['exnLen'] = sum(h['exnLenList'])
#
#     h['cdsList'] = []
#
#     for (s,e) in h['exnList']:
#
#         if s<=h['cdsSta'] and h['cdsSta']<=e:
#             s = h['cdsSta']
#
#         if s<=h['cdsEnd'] and h['cdsEnd']<=e:
#             e = h['cdsEnd']
#
#         if h['cdsSta']<=s and e<=h['cdsEnd']:
#             h['cdsList'].append((s,e))
#
#     h['cdsLen'] = sum([e-s for (s,e) in h['cdsList']])
#
#     return h
#
#
# def loadKgByChr(dataFileName='/Z/Sequence/ucsc_hg19/annot/knownGene.txt',h={}):
#
#     for line in open(dataFileName):
#
#         r = processKgLine(line)
#
#         jkbasic.addHash(h, r['chrom'], r)
#
#     return h
#
#
# def processKgLine(line):
#
#     tokL = line.rstrip().split('\t')
#
#     h = {}
#
#     h['geneId'] = tokL[0]
#     h['chrom'] = tokL[1]
#     h['chrNum'] = tokL[1][3:]
#     h['strand'] = tokL[2]
#     h['txnSta'] = int(tokL[3])
#     h['txnEnd'] = int(tokL[4])
#     h['txnLen'] = h['txnEnd'] - h['txnSta']
#     h['cdsSta'] = int(tokL[5])
#     h['cdsEnd'] = int(tokL[6])
#     h['exnList'] = map(lambda x,y: (int(x),int(y)), tokL[8].split(',')[:-1], tokL[9].split(',')[:-1])
#     h['exnLenList'] = [e-s for (s,e) in h['exnList']]
#     h['exnLen'] = sum(h['exnLenList'])
#
#     h['cdsList'] = []
#     frontL, backL = [],[]
#
#     for (s,e) in h['exnList']:
#
#         if e<=h['cdsSta']:
#             frontL.append((s,e))
#         elif s<=h['cdsSta'] and h['cdsSta']<=e:
#             frontL.append((s,h['cdsSta']))
#             h['cdsList'].append((h['cdsSta'],e))
#         elif h['cdsSta']<=s and e<=h['cdsEnd']:
#             h['cdsList'].append((s,e))
#         elif s<=h['cdsEnd'] and h['cdsEnd']<=e:
#             h['cdsList'].append((s,h['cdsEnd']))
#             backL.append((h['cdsEnd'],e))
#         elif h['cdsEnd']<=s:
#             backL.append((s,e))
#         else:
#             raise Exception
#
#     if h['strand'] == '+':
#         h['utr5pLen'] = sum([e-s for (s,e) in frontL])
#     elif h['strand'] == '-':
#         h['utr5pLen'] = sum([e-s for (s,e) in backL])
#     else:
#         raise Exception
#
#     h['cdsLen'] = sum([e-s for (s,e) in h['cdsList']])
#
#     if h['strand'] == '+':
#         exnLenListH = h['exnLenList']
#     else:
#         exnLenListH = h['exnLenList'][::-1]
#
#     transOffset = h['utr5pLen'] * -1
#     h['frame'] = []
#
#     # index of h['frame'] is [exon number]-1
#
#     for i in range(len(exnLenListH)):
#
#         if 0 <= transOffset < h['cdsLen']:
#             frame5p = transOffset % 3
#         else:
#             frame5p = None
#
#         if 0 <= transOffset+exnLenListH[i]-1 < h['cdsLen']:
#             frame3p = (transOffset+exnLenListH[i]-1) % 3
#         else:
#             frame3p = None
#
#         h['frame'].append((frame5p, frame3p))
#
#         transOffset += exnLenListH[i]
#
#     return h

def loadCosmic(cosmicDat='/data1/Sequence/cosmic/cosmic.dat'):
    h = {}

    for line in open(cosmicDat):
        colL = line.rstrip().split('\t')
        chr = colL[0]
        sta = colL[1]
        end = colL[2]
        ref = colL[4]
        alt = colL[5]
        key = (chr, sta, end, ref, alt)
        if key not in h:
            h[key] = 'Y'
    
    return h

def overlap(x,y): # s: base-0, e: base-1

    c1,s1,e1 = x
    c2,s2,e2 = y

    if c1 != c2 or e2<=s1 or e1<=s2:
        return 0

    s = max(s1,s2)
    e = min(e1,e2)

    if s < e:
        return e-s
    else:
        return 0


def mergeLoci(locusL,gap=10):

    if len(set([x.strand for x in locusL]))>1 or len(set([x.chrNum for x in locusL]))>1:
        print ('warning: heterogeneous mix of loci')
        print (locusL)
        return locusL

    locusL.sort(key=lambda x: x.chrEnd)
    locusL.sort(key=lambda x: x.chrSta)

    locusMergedL = []

    i = 0

    while i < len(locusL):

        chrS1, chrE1 = locusL[i].chrSta, locusL[i].chrEnd

        curE = chrE1

        j = i+1

        idL = [locusL[i].id]

        while j < len(locusL):

            chrS2, chrE2 = locusL[j].chrSta, locusL[j].chrEnd

            if curE + gap < chrS2:
                break

            curE = max(curE,chrE2)
            idL.append(locusL[j].id)

            j += 1

        newLocus = copy.deepcopy(locusL[i])
        newLocus.chrEnd = max(locusL[k].chrEnd for k in range(i,j))
        newLocus.id = '|'.join(idL)

        locusMergedL.append(newLocus)

        i = j

    return locusMergedL


def spliceAI_run(chrom,pos,ref,alt): # hg38, pos-1; indel exam: T to TA

    vcf = '''##fileformat=VCFv4.2
##fileDate=20191004
##reference=GRCh38/hg38
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr21,length=46709983>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
%s\t%s\t.\t%s\t%s\t.\t.\t.\n''' % (chrom,pos,ref,alt)

    import tempfile

    f = tempfile.NamedTemporaryFile()
    f.write(vcf); f.flush()

    return os.popen("cat %s | spliceai -R %s/D/Sequences/hg38/hg38.fa -A grch38" % (f.name,homedir)).readlines()



def variant_bi(chrNum,chrPos,strand,ref,alt,assembly='hg38',hexH4=False): # ref, alt: transcript base

    if not hexH4:
        hexH4 = load_hexamer4() 

    print ('Hexamer')
    print (hexamer4_byCoord_general(chrNum,chrPos-1,chrPos,strand,ref,alt,hexH4,assembly))
    print

    mes_result = mes_byCoord(chrNum,chrPos,strand,ref,alt,assembly,verbose=True) 

    print ('MES donor')
    print ('delta', map(float,list(mes_result[0])))
    print ('final', map(float,mes_result[2]))

    print ('MES acceptor')
    print ('delta', map(float,list(mes_result[1])))
    print ('final', map(float,mes_result[3]))

#def variant_ai(locusStr):
#
#    print spliceAI(locusStr)
#    print
#    print primateAI_raw(locusStr) 

def labranchor_query(locusStr):

    return tbi_bed_query('~/D/DB/labranchor/lstm.gencode_v19.hg38.top.sorted.bed.gz', locusStr)





#def hexamer_delta(s_before,s_after,hexH):
#
#    scores_after = hexamer(s_after,hexH)
#    scores_before = hexamer(s_before,hexH)
#
#    return scores_after - scores_before 

def hexamer_byCoord(chrNum,pos,strand,ref,alt,hexH,assembly,verbose=False):

    l = locus('chr%s:%s-%s%s' % (chrNum,pos-6,pos+5,strand))

    s_before = l.twoBitFrag(assembly)

    if s_before[5:5+len(ref)] != ref:
        print(s_before, ref, alt)
        raise Exception

    s_after = s_before[:5] + alt + s_before[-5:]

    scores_after = hexamer(s_after,hexH)
    scores_before = hexamer(s_before,hexH)

    if verbose:
        print(scores_before, scores_after, scores_after - scores_before)

    return sum(scores_after - scores_before)


def calColor(score):

    ceiling = 2

    if score>0:
        return (1,1-min(score,ceiling)/ceiling,1-min(score,ceiling)/ceiling)
    elif score<0:
        return (1-max(score,-ceiling)/-ceiling,1-max(score,-ceiling)/-ceiling,1)
    else:
        return (1,1,1)



class InitiationFailureException(Exception): pass



    
class transcript:


    def __init__(self,geneName,refFlatFileName='%s/Data/DB/refFlat_hg19.txt' % homedir,assembly='hg19'):  # return the longest transcript matching geneName

        rows = map(processRefFlatLine,os.popen('grep "^%s    " %s' % (geneName,refFlatFileName)).readlines())
        rows = filter(lambda x: not '_' in x['chrNum'], rows)

        if len(rows)==0:
            raise InitiationFailureException()

        rows.sort(lambda x,y: cmp(y['cdsLen'],x['cdsLen'])); row = rows[0] # take refSeq with longest coding region

        self.assembly = assembly

        self.geneName = row['geneName']
        self.refSeqId = row['refSeqId']
        self.chrom = row['chrom']
        self.chrNum = row['chrNum']
        self.strand = row['strand']

        self.txnLen = row['txnLen']

        self.exnList = row['exnList']
        self.exnLen = row['exnLen']

        self.cdsList = row['cdsList']
        self.cdsLen = row['cdsLen']

    def txnOverlap(self,region):

        return overlap((self.chrNum,self.exnList[0][0],self.exnList[-1][-1]),region)

    def cdsOverlap(self,region):

        total = 0

        for (cdsS,cdsE) in self.cdsList:

            total += overlap((self.chrNum,cdsS,cdsE),region)

        return total

    def exnOverlap(self,region):

        total = 0

        for (exnS,exnE) in self.exnList:

            total += overlap((self.chrNum,exnS,exnE),region)

        return total


def geneNameH(refFlatFileName='/Z/Sequence/ucsc_hg19/annot/refFlat.txt', knownToRefSeqFileName='/Z/Sequence/ucsc_hg19/annot/knownToRefSeq.txt', \
        hugoFileName='/Z/Sequence/geneinfo/hugo.txt'):

    geneNameH = {}

    for line in open(refFlatFileName):

        h = processRefFlatLine(line)

        geneNameH[h['refSeqId']] = h['geneName']
        geneNameH[h['geneName']] = h['geneName']

    for line in open(knownToRefSeqFileName):

        (knownId,refSeqId) = line[:-1].split('\t')

        try:
            geneNameH[knownId] = geneNameH[refSeqId]
        except:
            pass

    for line in open(hugoFileName):

        (geneName,geneDesc,aliases,geneCardNames,refSeqIds) = line[:-1].split('\t')

        for refSeqId in refSeqIds.split(','):
            
            if refSeqId not in geneNameH:
                geneNameH[refSeqId] = geneName

        for alias in aliases.split(','):

            if alias not in geneNameH:
                geneNameH[alias] = geneName

        for geneCardName in geneCardNames.split(','):

            geneNameH[geneCardName] = geneName

    return geneNameH


def geneSetH(biocartaFileName='/Z/Sequence/geneinfo/BIOCARTA.gmt', goFileName='/Z/Sequence/geneinfo/GO.gmt', keggFileName='/Z/Sequence/geneinfo/KEGG.gmt'):

    geneSetH = {'biocarta':{}, 'go':{}, 'kegg':{}}

    for line in open(biocartaFileName):

        tokL = line[:-1].split('\t')
        geneSetH['biocarta'][tokL[0]] = (tokL[1],tuple(tokL[2:]))

    for line in open(goFileName):

        tokL = line[:-1].split('\t')
        geneSetH['go'][tokL[0]] = (tokL[1],tuple(tokL[2:]))

    for line in open(keggFileName):

        tokL = line[:-1].split('\t')
        geneSetH['kegg'][tokL[0]] = (tokL[1],tuple(tokL[2:]))

    return geneSetH


def loadCensus(censusFileName='/Z/DB/cancer_gene_census_hg38.txt'):

    censusH = {}

    for line in open(censusFileName):

        tokL = line[:-1].split('\t')

        (geneName,desc,somatic,germline,role,mutType,translocPartners) = \
            (tokL[0],tokL[1],tokL[7],tokL[8],tokL[12],tokL[13],tokL[14])

        if geneName == 'Gene Symbol':
            continue

        censusH[geneName] = {'desc':desc}
        censusH[geneName]['somatic'] = somatic
        censusH[geneName]['germline'] = germline
        censusH[geneName]['role'] = role
        censusH[geneName]['mutType'] = mutType
        censusH[geneName]['translocPartners'] = translocPartners
    
    return censusH


def geneInfoH(geneNameH, geneSetH, refSeqSummaryFileName='/Z/Sequence/ucsc_hg19/annot/refSeqSummary.txt', hugoFileName='/Z/Sequence/geneinfo/hugo.txt', \
        censusFileName='/Z/Sequence/geneinfo/cancer_gene_census.txt', biocartaFileName='/Z/Sequence/geneinfo/BIOCARTA.gmt', \
        goFileName='/Z/Sequence/geneinfo/hugo.txt', keggFileName='/Z/Sequence/geneinfo/hugo.txt'):

    geneInfoH = {}

    for line in open(refSeqSummaryFileName):

        (refSeqId,status,summary) = line[:-1].split('\t')

        if refSeqId in geneNameH:

            geneName = geneNameH[refSeqId]

            if geneName not in geneInfoH:
                geneInfoH[geneName] = {}

            geneInfoH[geneName]['summary'] = summary

    for line in open(hugoFileName):

        (geneName,desc,aliases,geneCardName,refSeqIds) = line[:-1].split('\t')

        if geneName not in geneInfoH:
            geneInfoH[geneName] = {}

        geneInfoH[geneName]['desc'] = desc 
        geneInfoH[geneName]['aliases'] = aliases
        geneInfoH[geneName]['refSeqIds'] = refSeqIds

    for line in open(censusFileName):

        tokL = line[:-1].split('\t')

        (geneName,desc,somatic,germline,mutType,translocPartners) = (tokL[0],tokL[1],tokL[7],tokL[8],tokL[12],tokL[13])

        if geneName == 'Symbol':
            continue

        if geneName not in geneInfoH:
            geneInfoH[geneName] = {'desc':desc}

        geneInfoH[geneName]['census_somatic'] = somatic
        geneInfoH[geneName]['census_germline'] = germline
        geneInfoH[geneName]['census_mutType'] = mutType
        geneInfoH[geneName]['census_translocPartners'] = translocPartners


    for geneSetDB in geneSetH.keys():

        for (geneSetName,(geneSetDesc,geneNameL)) in geneSetH[geneSetDB].iteritems():

            for geneName in geneNameL:

                if geneName in geneInfoH:
                    jkbasic.addHash(geneInfoH[geneName],geneSetDB,(geneSetName,geneSetDesc))
                else:
                    geneInfoH[geneName] = {geneSetDB:[(geneSetName,geneSetDesc)]}

    return geneInfoH


class gene:


    def __init__(self,identifier,geneNameH=None,geneSetH=None,geneInfoH=None,geneDB={}):

        if geneDB != {}:

            if 'geneNameH' in geneDB:
                self.geneNameH = geneDB['geneNameH']
            else:
                self.geneNameH = geneNameH()

            if 'geneSetH' in geneDB:
                self.geneSetH = geneDB['geneSetH']
            else:
                self.geneSetH = geneSetH()

            if 'geneInfoH' in geneDB:
                self.geneInfoH = geneDB['geneInfoH']
            else:
                self.geneInfoH = geneInfoH(self.geneNameH,self.geneSetH)

            try:
                self.geneName = self.geneNameH[identifier]
            except:
                self.geneName = None

            if self.geneName and self.geneName in self.geneInfoH:
                self.geneInfo = self.geneInfoH[self.geneName]
            else:
                self.geneInfo = {}

        else:

            if geneNameH:
                self.geneNameH = geneNameH
            else:
                self.geneNameH = geneNameH()

            if geneSetH:
                self.geneSetH = geneSetH
            else:
                self.geneSetH = geneSetH()

            if geneInfoH:
                self.geneInfoH = geneInfoH
            else:
                self.geneInfoH = geneInfoH(geneNameH,geneSetH)

            try:
                self.geneName = self.geneNameH[identifier]
            except:
                self.geneName = None

            if self.geneName and self.geneName in geneInfoH:
                self.geneInfo = geneInfoH[self.geneName]
            else:
                self.geneInfo = {}

    def getAttr(self,attr):

        if attr in self.geneInfo:
            return self.geneInfo[attr]
        else:
            return ''

#def getGeneDB(geneNameH=geneNameH(), geneSetH=geneSetH(), geneInfoH=geneInfoH(geneNameH(),geneSetH())):
#    geneDB = {'geneNameH': geneNameH,'geneSetH': geneSetH, 'geneInfoH': geneInfoH}
#    return geneDB


def getFrameInfoH():

    kgH = loadKgByChr()
    frameInfoH = {}

    for chrom in kgH.keys():
        
        for t in kgH[chrom]:
            frameInfoH[t['geneId']] = t['frame']

    return frameInfoH


def frameCons(transId1,exnNum1,transId2,exnNum2,frameInfoH):
    
    if transId1 in frameInfoH:
        frame1 = frameInfoH[transId1][exnNum1-1][1]
    else:
        frame1 = None

    if transId2 in frameInfoH:
        frame2 = frameInfoH[transId2][exnNum2-1][0]
    else:
        frame2 = None

    if None not in (frame1,frame2):
        if ((2-frame1) + frame2) % 3 == 0:
            return 'Y'
        else:
            return 'N'
    else:
        return None

def lookupPileup(pileupDirL,sId,chrom,loc,ref,alt,flag='T'):

    inputFileNL = []

    if flag == 'T':
        for pileupDir in pileupDirL:
            inputFileNL += os.popen('find %s -name %s_T_*%s.pileup_proc' % (pileupDir,sId,chrom)).readlines()
    else:
        for pileupDir in pileupDirL:
            inputFileNL += os.popen('find %s -name %s_N_*%s.pileup_proc' % (pileupDir,sId,chrom)).readlines()
            inputFileNL += os.popen('find %s -name %s_B_*%s.pileup_proc' % (pileupDir,sId,chrom)).readlines()

    if len(inputFileNL) > 1:
        inputFileNL = filter(lambda x: not re.match('.*KN.*', x),inputFileNL)

    if len(inputFileNL) == 0:
        return None

    resultL = os.popen('grep -m 1 "^%s:%s," %s' % (chrom,loc,inputFileNL[0].rstrip()), 'r').readlines()

    if len(resultL)==0:
        return None
    else:
        tL = resultL[0].rstrip().split(',')
        if ref != tL[2]:
            sys.exit(1)
        refCount = int(tL[3])
        altCount = tL[4].count(alt)
        return (altCount,refCount)


## batch version of lookupPileup() without sample id
## output dictionary of 'sample':'refCount|altCount'
def lookupPileup_batch(pileupDirL,chrom,loc,ref,alt,flag='T',useFlag=True):

    ## same critera as mutScan
    minCover = 3
    minMutReads = 2
    minFreq = 0.01

    inputFileNL = []

    if useFlag:
        if flag == 'T':
            for pileupDir in pileupDirL:
                inputFileNL += os.popen('find %s -name *_T_*%s.pileup_proc' % (pileupDir,chrom)).readlines()
        else:
            for pileupDir in pileupDirL:
                inputFileNL += os.popen('find %s -name *_N_*%s.pileup_proc' % (pileupDir,chrom)).readlines()
                inputFileNL += os.popen('find %s -name *_B_*%s.pileup_proc' % (pileupDir,chrom)).readlines()
    else:
        for pileupDir in pileupDirL:
            inputFileNL += os.popen('find %s -name *%s.pileup_proc' % (pileupDir,chrom)).readlines()
    
    if len(inputFileNL) > 1:
        inputFileNL = filter(lambda x: not re.match('.*KN.*', x), inputFileNL)
    
    if len(inputFileNL) == 0:
        return None
    
    resultH = {}
    for inputFile in inputFileNL:
        sampN = inputFile.rstrip().split('/')[-1].split('_')[0]
        resultL = os.popen('grep -m 1 -P "^%s\\t%s\\t" %s' % (chrom, loc, inputFile.rstrip()), 'r').readlines()

        if len(resultL) > 0:
            tL = resultL[0].rstrip().split('\t')
            tot = int(tL[3])
            altCount = tL[4].upper().count(alt)
            refCount = tot - altCount
            resultH[sampN] = '%s|%s' % (refCount, altCount)
        else:
            resultH[sampN] = 'NA'

    return resultH

def df_gtex_tpm():
    
    # human gene expression
    
    return pd.read_csv('%s/D/DB/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz' % homedir,skiprows=2,sep='\t').iloc[:,1:].groupby('Description').median()

def df_exac():

    return pd.read_csv('%s/D/DB/exac/forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz' % homedir,sep='\t',dtype={'pLI':np.float}).iloc[:,1:].set_index('gene')

class tcgaCnaDB:

    def __init__(self,gctFileName):

        self.db= {}
        self.idx= {}

        inFile = open(gctFileName)

        inFile.readline(); inFile.readline()

        headerL = inFile.readline()[:-1].split('\t')

        for i in range(2,len(headerL)):
            self.idx[headerL[i]] = i-2

        for line in inFile:

            tokL = line[:-1].split('\t')
            self.db[tokL[0]] = tokL[2:]

    def query(self,sampN,geneN):

        if geneN in self.db and sampN in self.idx:
            return self.db[geneN][self.idx[sampN]]
        else:
            return ''



class locus: # UCSC type
    
    def __init__(self,loc,id=''):

        rm = re.match('([^:]+):([0-9,]+)-([0-9,]+)([+-])',loc) # base-0, base-1

        if rm:

            self.strand = rm.group(4)
            self.chrom = rm.group(1)

            if self.chrom[:3] == 'chr':
                self.chrNum = rm.group(1)[3:]
            else:
                self.chrNum = None

            self.chrSta = int(rm.group(2))
            self.chrEnd = int(rm.group(3))

            self.id = id

        else:

            rm = re.match('([+-])([^:]+):([0-9,]+)..([0-9,]+)',loc) # base-1, base-1

            if rm:

                self.strand = rm.group(1)
                self.chrom = rm.group(2)

                if self.chrom[:3] == 'chr':
                    self.chrNum = rm.group(2)[3:]
                else:
                    self.chrNum = None

                chrPosL = [int(rm.group(3)), int(rm.group(4))]

                self.chrSta = min(chrPosL) - 1
                self.chrEnd = max(chrPosL) 

                self.id = id

            else:

                raise Exception

    def toString(self,style='UCSC'):

        if style=='gsnap':
            return '%s%s:%s..%s' % (self.strand,self.chrom,self.chrSta+1,self.chrEnd)
        else:
            return '%s:%s-%s%s' % (self.chrom,self.chrSta,self.chrEnd,self.strand)

    def overlap(self,region):

        return overlap((self.chrom,self.chrSta,self.chrEnd),region)

    def overlappingGeneL(self,refFlatH=None,refFlatFileName='~/D/Sequences/ucsc_hg19_ref/refseq',strand_sensitive=False):

        gL = set()

        if refFlatH == None and refFlatFileName != '':
            refFlatH = loadRefFlatByChr(refFlatFileName)
        
        if self.chrom not in refFlatH:
            return []

        for l in refFlatH[self.chrom]:

            if strand_sensitive:

                if self.overlap((l['chrom'],l['txnSta'],l['txnEnd'])) > 0 and self.strand==l['strand']:
                    gL.add(l['geneName'])

            else:

                if self.overlap((l['chrom'],l['txnSta'],l['txnEnd'])) > 0:
                    gL.add(l['geneName'])

        return tuple(gL)

    def regionType(self,h=loadBlatOutputByChr()): 
        return getRegionType(h,self)


    def twoBitFrag(self, assembly='hg38', buffer5p=0, buffer3p=0):
        import platform as pl
        sep = os.path.sep
        homedir = os.path.abspath(os.path.dirname(__file__)+sep+'data')
        twoBitFilePath='%s/D/Sequences/%s/%s.2bit' % (homedir,assembly,assembly)

        if self.strand == '+':
            staPos = self.chrSta - buffer5p
            endPos = self.chrEnd + buffer3p
        else:
            staPos = self.chrSta - buffer3p
            endPos = self.chrEnd + buffer5p
        print('  chrom:',self.chrom,'  stapos: ',staPos,'  endpos: ',endPos)
        if pl.system().lower() =='windows':
            fragFile = os.popen
        else:            
            fragFile = os.popen('%s/tools/jkent/twoBitToFa %s:%s:%s-%s stdout' % (homedir, twoBitFilePath, self.chrom, staPos, endPos), 'r')
        fragFile.readline()

        seq = fragFile.read().replace('\n','').rstrip().upper()

        if self.strand == '+':
            return seq
        else:
            return jkbio.rc(seq)
def getRegionType(h,loc): # locusTupe = (chrom,sta,end) 1-base

    locT = (loc.chrom,loc.chrSta+1,loc.chrEnd)
    regions = []
    for t in h[locT[0]]:

        if t['txnEnd'] < locT[1]:
            continue

        elif locT[2] <= t['txnSta']:
            break

        if t['strand']==loc.strand:
            sense = 'sense'
        else:
            sense = 'antisense'

        flag = None

        for regionName in ['cdsList','utr5','utr3','intron']:

            for s,e in t[regionName]:

                if overlap(('_',locT[1]-1,locT[2]),('_',s,e)) == 0:
                    continue

                if regionName == 'cdsList':

                    if locT[2]-locT[1] == 0:

                        offset = 0

                        for s,e in t['cdsList']:
                            offset += min(e,locT[1])-min(s,locT[1])

                        frame = (offset-1) % 3

                        if t['strand'] == '-':
                            frame= 2-frame

                        flag = 'cds_%d_%d' % (offset-1,frame)

                    else:

                        flag = 'cds_m'

                else:

                    flag = regionName

                marg = margin((locT[1]-1,locT[2]),(s,e))

                if (marg[0] or marg[1]) < 0:
                    continue

                else:
                    regions.append((t['transName'],t['transID'],flag, sense, marg))
        if t['cdsList'] == []:
            regions.append((t['transName'],t['transID'],'lnc', sense, -1))

    return list(set(regions))