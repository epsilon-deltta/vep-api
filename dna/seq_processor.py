    
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
# =========================================================
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
    print(l_sign,alt_sign,r_sign)
    # ref > alt
    # alt < ref
    new_context = []
    print("alt_sign",alt_sign)
    for chr_info in context:

        if '<'   in alt_sign :
            if l_sign == chr_info['alt'] and r_sign == chr_info['ref']:
                new_context.append(chr_info)
        elif '>' in alt_sign:
            print("elif")
            print("r  :",r_sign,"   l:",l_sign)
            print("ref :",chr_info['ref'])
            if l_sign == chr_info['ref'] and r_sign == chr_info['alt']:
                print("if")
                new_context.append(chr_info)

    return new_context 
