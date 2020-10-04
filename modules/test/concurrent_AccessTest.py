import sys ,os
sys.path.append(os.path.dirname(__file__)+'/..')
import jkgenome as jk
import threading as thr
import time
#ca (concurrent accesses)
class caTest():
    def __init__(self):
        self.locus_str    = '11:108236168-10823616'
        self.bedFilePath  = '%s/BigFiles/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz' % jk.homedir
    def test_tbi_bed_query(self):
        ths=[]
        for _ in range(10):
            ths.append( thr.Thread(target=jk.tbi_bed_query, args=(self.bedFilePath,self.locus_str,) ) )
        for n in range(10):
            ths[n].start() 
                    
if __name__ == "__main__":
    test = caTest()
    test.test_tbi_bed_query()
    