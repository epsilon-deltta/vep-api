from unittest import TestCase, main
from jkgenomeTest import jkgTest as jt
class jkg_trackTest(TestCase):    
    
    def test_spliceAI_track(self,seq='11:108236168-10823616'):
        jt.test_spliceAI()
        jt.test_spliceAI_raw()
        jt.test_tbi_bed_query()

if __name__ == '__main__':
    main()