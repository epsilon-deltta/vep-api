    
import base64
import os

import jkgenome as jk

def get_mes3(seq):
    return jk.mes3(seq)
# base64    
def get_hexamer_track_info(seq='GCACCAAAGGAGGAAATCTGCCTGCATGATCCAATCACCTCCCATCAGGCCCCACCTCCAACATTGGGTATT'):
    seq.upper()
    if seq.count(',') == 3:
        seq = get_position_To_seq(seq)

    # print("="*20,'\n',"="*20,'\n',"="*20,'\n',"="*20)
    # print("this is seq :",seq)
    track_data        = jk.hexamer_track(seq)
    track_data['seq'] = seq
    fpath             = track_data['path']
    
    if os.path.isfile(fpath):
        with open(fpath, 'rb') as img:
            track_data['base64'] = base64.b64encode(img.read() )
        os.remove(fpath)
        print("deleted file")
    else:
        print("File not Found!!!!")
    del track_data['path']
    return track_data

def get_position_To_seq(seq="chr14,31549779,31549850,+"):
    arg   = seq.split(',')
    locus = jk.locus('%s:%s-%s%s' % tuple(arg))
    seq   = locus.twoBitFrag('hg38')
    return seq
# =================splice AI========================================
# input : 11:108236168-108236168 A>G 
#output :[ {'chrN': '11', 'pos': 108236168, 'ref': 'AGTAG', 'alt': 'A', 'geneName': 'ATM', 'score': {'AG': 0.0, 'AL': 0.0, 'DG': 0.02, 'DL': 0.09}, 'relPos': {'AG': 49, 'AL': 0, 'DG': -48, 'DL': 0}},..]
def get_splice_ai(seq='11:108236168-108236168'): 
    print("seq str is ",seq)
    seq = seq.rstrip()
    filter_sign = None 


    if ' ' in seq:
        filter_sign = seq.split(' ')[1]
        seq         = seq.split(' ')[0]

    if '-' not in seq:
        seq = onePos_To_twoPos(seq)
    context = jk.spliceAI(seq)

    if filter_sign is not None :
        context = get_matched_seq(context,filter_sign) 
    return context
def onePos_To_twoPos(seq='11:108236168'):
    second_pos = seq.split(':')[1]
    seq        = seq + '-' + second_pos
    return seq
import re
def get_matched_seq(context,filter_sign):
    signs = re.split(r"[/>/<]",filter_sign)

    l_sign   = signs[0].upper()
    alt_sign = filter_sign[len(l_sign)]
    r_sign   = signs[1].upper()
    # print(l_sign,alt_sign,r_sign) for debuging
    # ref > alt
    # alt < ref
    new_context = []
    for chr_info in context: 

        if '<'   in alt_sign :
            if l_sign == chr_info['alt'] and r_sign == chr_info['ref']:
                new_context.append(chr_info)
        elif '>' in alt_sign:
            if l_sign == chr_info['ref'] and r_sign == chr_info['alt']:
                new_context.append(chr_info)

    return new_context 
# =================splice AI pipeline========================
# variant_bi2
def get_splice_ai_option(seq,options="hex_mes"):
    print("pass")
    splice_list  = get_splice_ai(seq)
    print(splice_list)
    handled_list = handle_option(splice_list,options)
    return handled_list
def handle_option(splice_list,options):
    if options =="hex_mes":
        for i, item in enumerate(splice_list):
            hex_mes = get_hexamer_maxent(item)
            splice_list[i].update(hex_mes)
    return splice_list
# input :{'chrN': '11', 'pos': 108236168, 'ref': 'AGTAG', 'alt': 'A', 'geneName': 'ATM', 'score': {'AG': 0.0, 'AL': 0.0, 'DG': 0.02, 'DL': 0.09}, 'relPos': {'AG': 49, 'AL': 0, 'DG': -48, 'DL': 0}}
# output :
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
# =========================================
    
def get_hexamer_maxent(sn):
    return jk.variant_bi2(sn['chrN'],sn['pos'],sn['ref'],sn['alt'])
# input  11:108236168 A>C
# output :output.keys() : chrN, pos,ref,alt , (variant_bi2 dict.keys())
def get_hexamer_maxent1(seq):
    seq = seq.rstrip()

    sn               = {}
    seq ,filter_sign = seq.split(' ')

    sn['chrN']  = int(seq.split(':')[0] )
    sn['pos']   = int(seq.split(':')[1] )

    signs = re.split(r"[/>/<]",filter_sign)

    l_sign   = signs[0].upper()
    alt_sign = filter_sign[len(l_sign)]
    r_sign   = signs[1].upper()

    if '<'  in alt_sign :
        sn['alt']   = l_sign 
        sn['ref']   = r_sign
    elif '>' in alt_sign :
        sn['alt']   = r_sign
        sn['ref']   = l_sign
    print("sn is",sn)
    sn.update( get_hexamer_maxent(sn) )
    return sn

