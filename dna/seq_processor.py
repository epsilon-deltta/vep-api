import jkgenome as jk
import base64
import os
def get_mes3(seq):
    return jk.mes3(seq)
# base64    
def get_hexamer_track_info(seq='GCACCAAAGGAGGAAATCTGCCTGCATGATCCAATCACCTCCCATCAGGCCCCACCTCCAACATTGGGTATT'):
    seq.upper()
    track_data = jk.hexamer_track(seq)
    fpath = track_data['path']
    
    if os.path.isfile(fpath):
        with open(fpath, 'rb') as img:
            track_data['base64'] = base64.b64encode(img.read() )
        os.remove(fpath)
        print("deleted file")
    else:
        print("File not Found!!!!")
    del track_data['path']
    return track_data
