import sys ,os
sys.path.append(os.path.dirname(__file__)+'/..')
import jkgenome as jk
import unittest
homedir = jk.homedir
class jkgTest(unittest.TestCase):    

    def test_spliceAI(self,seq='11:108236168-10823616'):
        print("=====spliceAI Test")
        print("func : spliceAI")
        result = jk.spliceAI('11:108236168-10823616')
        if result == []:
            print("splice ai result is empty")
            print("result :",result)
        
    def test_spliceAI_raw(self,seq='11:108236168-10823616'):
        print("func : spliceAI_raw")
        result = jk.spliceAI_raw(seq)
        
        if result == []:
            print("spliceAI_raw's result is empty")
            print("result :",result)
    def test_tbi_bed_query(self,bedFilePath='%s/BigFiles/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz' % homedir
                               ,locusStr='11:108236168-10823616'):
        print("func : tbi_bed_query")
        result = jk.tbi_bed_query(bedFilePath,locusStr)
        
        if result == []:
            print("spliceAI_raw's result is empty")
            print("result :",result)
            
            # print("inner cmd cheking : ")
    def test_spliceAI_track(self,seq='11:108236168-10823616'):
        self.test_spliceAI()
        self.test_spliceAI_raw()
        self.test_tbi_bed_query()

if __name__ == '__main__':
    jkt = jkgTest()
    # jkt.test_spliceAI()
    # jkt.test_spliceAI_raw()
    # jkt.test_tbi_bed_query()
    jkt.test_spliceAI_track()
    
    
        

        
