# fig.savefig('%shexamer_track.pdf' % (os.path.dirname(__file__)+sp+'data'+sp) )
import jkgenome as jk
seq ='GCACCAAAGGAGGAAATCTGCCTGCATGATCCAATCACCTCCCATCAGGCCCCACCTCCAACATTGGGTATT'
print(jk.hexamer_track(seq))