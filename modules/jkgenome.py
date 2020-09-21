
import random as rd
def mes3(seq): # seq example: ttcca aacga acttt tgtag GGA (23)
    try:
        if len(seq) != 23:
            raise Exception
    except Exception:
        pass
    # return float(os.popen('cd %s/tools/maxentscan; perl score3_mod.pl %s' % (homedir,seq)).readline())
    return score3_mod(seq)

#####=method
#####stub object 
def score3_mod(seq):
    return rd.random()
