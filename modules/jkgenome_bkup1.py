# mygenome

import sys, os, copy, re, gzip, collections
import jkbasic, jkbio


#RTK = ['AATK','AATYK','AATYK2','AATYK3','ACH','ALK','anaplastic lymphoma kinase','ARK','ATP:protein-tyrosine O-phosphotransferase [ambiguous]','AXL','Bek','Bfgfr','BRT','Bsk','C-FMS','CAK','CCK4','CD115','CD135','CDw135','Cek1','Cek10','Cek11','Cek2','Cek3','Cek5','Cek6','Cek7','CFD1','CKIT','CSF1R','DAlk','DDR1','DDR2','Dek','DKFZp434C1418','Drosophila Eph kinase','DRT','DTK','Ebk','ECK','EDDR1','Eek','EGFR','Ehk2','Ehk3','Elk','EPH','EPHA1','EPHA2','EPHA6','EPHA7','EPHA8','EPHB1','EPHB2','EPHB3','EPHB4','EphB5','ephrin-B3 receptor tyrosine kinase','EPHT','EPHT2','EPHT3','EPHX','ERBB','ERBB1','ERBB2','ERBB3','ERBB4','ERK','Eyk','FGFR1','FGFR2','FGFR3','FGFR4','FLG','FLK1','FLK2','FLT1','FLT2','FLT3','FLT4','FMS','Fv2','HBGFR','HEK11','HEK2','HEK3','HEK5','HEK6','HEP','HER2','HER3','HER4','HGFR','HSCR1','HTK','IGF1R','INSR','INSRR','insulin receptor protein-tyrosine kinase','IR','IRR','JTK12','JTK13','JTK14','JWS','K-SAM','KDR','KGFR','KIA0641','KIAA1079','KIAA1459','Kil','Kin15','Kin16','KIT','KLG','LTK','MCF3','Mdk1','Mdk2','Mdk5','MEhk1','MEN2A/B','Mep','MER','MERTK','MET','Mlk1','Mlk2','Mrk','MST1R','MTC1','MUSK','Myk1','N-SAM','NEP','NET','Neu','neurite outgrowth regulating kinase','NGL','NOK','nork','novel oncogene with kinase-domain','Nsk2','NTRK1','NTRK2','NTRK3','NTRK4','NTRKR1','NTRKR2','NTRKR3','Nuk','NYK','PCL','PDGFR','PDGFRA','PDGFRB','PHB6','protein-tyrosine kinase [ambiguous]','protein tyrosine kinase [ambiguous]','PTK','PTK3','PTK7','receptor protein tyrosine kinase','RET','RON','ROR1','ROR2','ROS1','RSE','RTK','RYK','SEA','Sek2','Sek3','Sek4','Sfr','SKY','STK','STK1','TEK','TIE','TIE1','TIE2','TIF','TKT','TRK','TRKA','TRKB','TRKC','TRKE','TYK1','TYRO10','Tyro11','TYRO3','Tyro5','Tyro6','TYRO7','UFO','VEGFR1','VEGFR2','VEGFR3','Vik','YK1','Yrk']
#
#TK = RTK + ['ABL','ABL1','ABL2','ABLL','ACK1','ACK2','AGMX1','ARG','ATK','ATP:protein-tyrosine O-phosphotransferase [ambiguous]','BLK','Bmk','BMX','BRK','Bsk','BTK','BTKL','CAKb','Cdgip','CHK','CSK','CTK','CYL','cytoplasmic protein tyrosine kinase','EMT','ETK','Fadk','FAK','FAK2','FER','Fert1/2','FES','FGR','focal adhesion kinase','FPS','FRK','FYN','HCK','HCTK','HYL','IMD1','ITK','IYK','JAK1','JAK2','JAK3','Janus kinase 1','Janus kinase 2','Janus kinase 3','JTK1','JTK9','L-JAK','LCK','LSK','LYN','MATK','Ntk','p60c-src protein tyrosine kinase','PKB','protein-tyrosine kinase [ambiguous]','PSCTK','PSCTK1','PSCTK2','PSCTK4','PSCTK5','PTK2','PTK2B','PTK6','PYK2','RAFTK','RAK','Rlk','Sik','SLK','SRC','SRC2','SRK','SRM','SRMS','STD','SYK','SYN','Tck','TEC','TNK1','Tsk','TXK','TYK2','TYK3','YES1','YK2','ZAP70']

assemblyH = {'hg18':'/data1/Sequence/ucsc_hg18/hg18_nh.fa', 'hg19':'/data1/Sequence/ucsc_hg19/hg19_nh.fa', 'hg38':'/Z/DB/Homo_sapiens.GRCh38.dna.primary_assembly.fa'}

import socket

if socket.gethostname() in ['cpu1','cpu2','gpu1','gpu2']:
    homedir = '/home/jinkuk/'
else:
    homedir = '/Users/jinkuk/'


def loadFasta(fastaPath): # only when seqs are in the same line

    h = collections.defaultdict(list)

    if fastaPath.endswith('.gz'):
        f = gzip.open(fastaPath)
    else:
        f = open(fastaPath)

    while 1:

        header = f.readline()
        
        if not header:
            break

        if not header.startswith('>'):
            raise Exception
            
        h[header.split('|')[0][1:]].append(f.readline().rstrip())

    return h


def loadBlatOutput(blatOutputPath='%s/D/Sequences/gencode_24/wgEncodeGencodeCompV24lift37.txt.gz' % homedir,by='transID'):

    h = collections.defaultdict(list)

    if blatOutputPath.endswith('.gz'):
        f = gzip.open(blatOutputPath)
    else:
        f = open(blatOutputPath)

    for line in f:
    
        r = processBlatLine(line)
        h[r[by]].append(r)

    from operator import itemgetter, attrgetter

    for k,vL in h.iteritems():
        h[k] = sorted(vL,key=itemgetter('txnSta','txnEnd'))

    return h

def loadBlatOutputByChr(blatOutputPath):

    return loadBlatOutput(blatOutputPath, 'chrom')

def loadBlatOutputByID(blatOutputPath):

    return loadBlatOutput(blatOutputPath,'transID')

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
    h['exnList'] = map(lambda x,y: (int(x),int(y)), tokL[9].split(',')[:-1], tokL[10].split(',')[:-1])
    h['exnLenList'] = [e-s for (s,e) in h['exnList']]
    h['exnLen'] = sum(h['exnLenList'])

    h['cdsList'] = []
    frontL, backL = [],[]

    if h['cdsSta'] != h['cdsEnd']:

        for (s,e) in h['exnList']:

            frontL.append((min(s,h['cdsSta']),min(e,h['cdsSta'])))
            h['cdsList'].append((max(s,h['cdsSta']),min(e,h['cdsEnd'])))
            backL.append((max(s, h['cdsEnd']), max(e, h['cdsEnd'])))

        frontL = filter(lambda x: x[0] < x[1],frontL)
        h['cdsList'] = filter(lambda x: x[0] < x[1], h['cdsList'])
        backL = filter(lambda x: x[0] < x[1], backL)

    if h['strand'] == '+':
            h['utr5'] = frontL
            h['utr3'] = backL
    elif h['strand'] == '-':
            h['utr5'] = backL
            h['utr3'] = frontL
    else:
            raise Exception

    return h

def getRegionType(h,locT): # locusTupe = (chrNum,sta,end) 1-base

    regions = []

    for t in h[locT[0]]:

        if t['txnEnd'] < locT[1]:
            continue

        elif locT[2] <= t['txnSta']:
            break

        flag = ''

        for regionName in ['cdsList','utr5','utr3']:
            if sum([overlap(('_',locT[1]-1,locT[2]),('_',s,e)) for s,e in t[regionName]]) > 0:
                flag = regionName
                break

        if flag == 'cdsList':

            if locT[2]-locT[1] == 0:

                offset = 0

                for s,e in t['cdsList']:
                    offset += min(e,locT[1])-min(s,locT[1])

                offset = (offset-1) % 3

                if t['strand'] == '-':
                    offset = 2-offset

                flag = 'cds_%d' % offset

            else:
                flag = 'cds_m'

        elif flag == '':

            if sum([overlap(('_',locT[1]-1,locT[2]),('_',s,e)) for s,e in t['exnList']]) > 0:
                flag = 'lnc'
            else:
                flag = 'intr'

        regions.append(flag)

    return list(set(regions))

def convertTrans2Genome(blatH,transID,transPos,transLen): # pos base-1

    result = []

    for trans in blatH[transID]:

        if '_' in  trans['chrom']:
            continue

        if transLen != trans['exnLen']:
            print(transLen, trans['exnLen'])
            raise Exception

        if trans['strand'] == '-':
            transPos = transLen-transPos+1

        tally = 1

        for i,n in enumerate(trans['exnLenList']):

            if tally <= transPos < tally+n:
                result.append((trans['chrom'], trans['exnList'][i][0]+(transPos-tally)+1, trans['strand'])) # pos base-1
                break

            tally += n

    return result

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

    (c1,s1,e1),(c2,s2,e2) = x,y

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
        print( 'warning: heterogeneous mix of loci', locusL)
        return locusL

    locusL.sort(lambda a,b: cmp(a.chrEnd,b.chrEnd))
    locusL.sort(lambda a,b: cmp(a.chrSta,b.chrSta))

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


class InitiationFailureException(Exception): pass


class locus: # UCSC type

    def __init__(self,loc,id=''):

        rm = re.match('([^:]+):([0-9,]+)-([0-9,]+)([+-])',loc)

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

            rm = re.match('([+-])([^:]+):([0-9,]+)..([0-9,]+)',loc)

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

    def overlappingGeneL(self,refFlatH=None,refFlatFileName='/Z/Sequence/ucsc_hg19/annot/refFlat.txt',strand_sensitive=False):

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

    def nibFrag(self, nibFragBase='/%s/D/Sequences/hg19/refFlat_hg19.txt' % (homedir,), buffer5p=0, buffer3p=0):

        if self.strand == '+':
            staPos = self.chrSta - buffer5p
            endPos = self.chrEnd + buffer3p
        else:
            staPos = self.chrSta - buffer3p
            endPos = self.chrEnd + buffer5p

        nibFragFile = os.popen('nibFrag -name="" %s/%s.nib %s %s %s stdout' % (nibFragBase, self.chrom, staPos, endPos, self.strand), 'r')

        nibFragFile.readline()

        return nibFragFile.read().replace('\n','').rstrip().upper()

    def twoBitFrag(self, twoBitFilePath='/%s/D/Sequences/hg19/hg19.2bit' % (homedir,), buffer5p=0, buffer3p=0):

        if self.strand == '+':
            staPos = self.chrSta - buffer5p
            endPos = self.chrEnd + buffer3p
        else:
            staPos = self.chrSta - buffer3p
            endPos = self.chrEnd + buffer5p

        fragFile = os.popen('%s/Z/tools/kentutils/twoBitToFa %s:%s:%s-%s stdout' % (homedir, twoBitFilePath, self.chrom, staPos, endPos), 'r')
        fragFile.readline()

        seq = fragFile.read().replace('\n','').rstrip().upper()

        if self.strand == '+':
            return seq
        else:
            return jkbio.rc(seq)

class transcript:


    def __init__(self,geneName,refFlatFileName='/Users/jinkuk/Data/DB/refFlat_hg18.txt',assembly='hg18'):  # return the longest transcript matching geneName

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
