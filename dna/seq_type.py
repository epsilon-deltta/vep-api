import re
import jkgenome as jk
# input > output
# input : chr1:1234567890 >False
# input : 1:1234567890 >False
# input : NM_000552.4(VWF):c.3797C>A >True
# input : NM_000552:3797C>A >True
# input : NM_000552:3797
def check_transID_Type(seq):
    # chr1:1234567890 >False
    if seq[:3] == 'chr' :
        return False
    # 1:1234567890 >False
    if bool(re.match(r"\d",seq ) ):
        return False
    return True
# input : NM_000552.4(VWF):c.3797C>A 
# input : NM_000552:3797C>A 
# input : NM_000552:3797 
# output : 12:6022032 C>A
def get_transId2pos(seq):
    transID  =  seq.split(':')[0]

    if '.' in transID :
        transID = transID[:transID.find('.')]
    transPos =  re.search(r"(?P<pos>[a-zA-Z]+[.]\d+|\d+)" , seq.split(':')[1] )
    transPos =  transPos.group('pos')

    if '.' in transPos :
        transPos = transPos[transPos.find('.')+1:]

    variant  =  seq.split(':')[1][len(transPos):]
    if exist_Variant(variant):
        variant = extract_variant(variant)
    else :
        variant = ''
    # output :[('chr12', 6022032, '-')]
    chrN ,pos, _ =  jk.convertTrans2Genome(transID, int(transPos) )[0]

    position = chrN+":"+str(pos)+' '+variant
    
    return position.rstrip()

def extract_variant(seq):
    result = re.search(r'[a-zA-Z]+[\>\<]+[a-zA-Z]+',seq)
    if result is  None:
        return ''
    return result.group(0)
def exist_Variant(seq):
    return bool(extract_variant(seq) )


# check_transID_Type("0ca2134")