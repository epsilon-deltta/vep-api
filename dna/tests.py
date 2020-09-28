from django.test import TestCase

# Create your tests here.
import sys, os
sep = os.path.sep
module_path = os.path.dirname(__file__)+sep+".."+sep+"modules"+sep
print("module_path :",module_path)
sys.path.append(module_path)
import seq_processor as sp
def unitTest():
    mes3_str = 'TTTCTGTGTCCTTCCTCCAACCC'
    print("==get_mes3==")
    print("input  : ",mes3_str)
    # print("output : ",sp.get_mes3(mes3_str))
    print("="*10)

    test_str ="GCACCAAAGGAGGAAATCTGCCTGCATGATCCAATCACCTCCCATCAGGCCCCACCTCCAACATTGGGTATT"
    print("==get_hexamer_track_info==")
    print("input",test_str)
    print("output ",sp.get_hexamer_track_info(test_str)['5p_exon'])
    print("="*10)


    region_str ='chr14,31549779,31549850,+'
    print("==get_position_to_seq==")
    print("input ",region_str)
    print("output ",sp.get_position_To_seq(region_str))
    print("="*10)

unitTest()