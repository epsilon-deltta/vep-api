import sys ,os
sys.path.append(os.path.dirname(__file__)+'/..')
import jkgenome as jk
from colorama import Fore ,Style
import unittest
homedir = jk.homedir
class jkgTest(unittest.TestCase):    
    # def setUp(self):
    #     # 테스트 시작 전 셋업한다. 각 test마다 매번 호출된다
    #     print(f'{Fore.YELLOW}================SET UP================ {Style.RESET_ALL}')
    def test_spliceAI(self,seq='11:108236168-10823616'):
        print("=====spliceAI Test")
        print(f"{Fore.GREEN}func : spliceAI{Style.RESET_ALL}")
        result = jk.spliceAI('11:108236168-10823616')
        if result == []:
            print("splice ai result is empty")
            print("result :",result)
        
    def test_spliceAI_raw(self,seq='11:108236168-10823616'):
        print(f"{Fore.GREEN}func : spliceAI_raw{Style.RESET_ALL}")
        result = jk.spliceAI_raw(seq)
        
        if result == []:
            print("spliceAI_raw's result is empty")
            print("result :",result)
    def test_tbi_bed_query(self,bedFilePath='%s/BigFiles/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz' % homedir
                               ,locusStr='11:108236168-10823616'):
        print(f"{Fore.GREEN}func : tbi_bed_query{Style.RESET_ALL}")
        result = jk.tbi_bed_query(bedFilePath,locusStr)
        
        if result == []:
            print("spliceAI_raw's result is empty")
            print("result :",result)
            
            # print("inner cmd cheking : ")

    # def tearDown(self):
    #     # 테스트 완료 후 매번 실행 
    #     print(f'{Fore.YELLOW}================Tear down================ {Style.RESET_ALL}')
if __name__ == '__main__':
    jkt = jkgTest()
    # jkt.test_spliceAI()
    # jkt.test_spliceAI_raw()
    # jkt.test_tbi_bed_query()
    # jkt.test_spliceAI_track()
    unittest.main()
    
    
        

        
