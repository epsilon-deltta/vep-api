# fig.savefig('%shexamer_track.pdf' % (os.path.dirname(__file__)+sp+'data'+sp) )
import jkgenome 
# seq ='GCACCAAAGGAGGAAATCTGCCTGCATGATCCAATCACCTCCCATCAGGCCCCACCTCCAACATTGGGTATT'
sd="sdfsa,dfas,dfs,fsdf,"
if sd.count(',') ==4:
    print("yes this is 4")
arg = seq.split(',')
locus = jkgenome.locus('%s:%s-%s%s' % tuple(arg))
seq   = locus.twoBitFrag('hg38')