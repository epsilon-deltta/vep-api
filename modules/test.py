# fig.savefig('%shexamer_track.pdf' % (os.path.dirname(__file__)+sp+'data'+sp) )
import jkgenome 
# seq ='GCACCAAAGGAGGAAATCTGCCTGCATGATCCAATCACCTCCCATCAGGCCCCACCTCCAACATTGGGTATT'
# print(jk.hexamer_track(seq))
print("testing")
region = ('chr14',31549779,31549850,'+')
print("testing")
locus = jkgenome.locus('%s:%s-%s%s' % region) #?
print("testing")
seq = locus.twoBitFrag('hg38')
print("testing")
len(seq)
print(seq)