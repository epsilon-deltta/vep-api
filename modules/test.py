import jkgenome as jk
import os
import time
print("=====start")
# result = jk.spliceAI('11:108241292-108241292') #output 01
def spliceai_test():    
    result = jk.spliceAI('11:108236168-108236168') #output 01
    for x in result:
        print(x)
def variant_bi2_test():
    jk.variant_bi2(12,51391349,'T','G')
if __name__ == '__main__':
    variant_bi2_test()
# /HDD8/ypil/djg-api-test/modules/data/tools/tabix-0.2.6/tabix /HDD8/ypil/djg-api-test/modules/data/BigFiles/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz 11:108236168-10823616
