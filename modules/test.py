import jkgenome as jk
import os
import time
print("=====start")
# result = jk.spliceAI('11:108241292-108241292') #output 01
result = jk.spliceAI('11:108236168-108236168') #output 01

# print(result)

for x in result:
    print(x)
# /HDD8/ypil/djg-api-test/modules/data/tools/tabix-0.2.6/tabix /HDD8/ypil/djg-api-test/modules/data/BigFiles/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz 11:108236168-10823616
