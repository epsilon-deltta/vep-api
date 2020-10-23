import vep_request as vr
def get_vep_data_test(seq='11:108236168 A>C'):
    vr.get_vep_data(seq)

if __name__ == '__main__':
    get_vep_data_test()