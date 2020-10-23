import os
# input 11:1234455 A>G
# output [{ ... } , { ... }]
# [{
#     "Uploaded_variation": "11_108236168_A/C",
#     "Location": "11:108236168",
#     "Allele": "C",
#     "Gene": "ENSG00000149311",
#     "Feature": "ENST00000278616",
#     "Feature_type": "Transcript",
#     "Consequence": "intron_variant",
#     "cDNA_position": "-",
#     "CDS_position": "-",
#     "Protein_position": "-",
#     "Amino_acids": "-",
#     "Codons": "-",
#     "Existing_variation": "-",
#     "Extra": "IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=SNV;SYMBOL=ATM;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:795;BIOTYPE=protein_coding;CANONICAL=YES;TSL=5;APPRIS=P1;CCDS=CCDS31669.1;ENSP=ENSP00000278616;SWISSPROT=Q13315;TREMBL=A0A024R3C7;UNIPARC=UPI000016B511;GENE_PHENO=1;INTRON=5/62\n"
# },...]
def get_vep_data(seq,option ='-e'):
    chrom , pos,ref,alt = split_seq(seq)
    vcf = '''##fileformat=VCFv4.2
##fileDate=20191004
##reference=GRCh38/hg38
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr21,length=46709983>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
%s\t%s\t.\t%s\t%s\t.\t.\t.\n''' % (chrom,pos,ref,alt)
    print('%s\t%s\t.\t%s\t%s\t.\t.\t.\n' % (chrom,pos,ref,alt))
    import tempfile

    f = tempfile.NamedTemporaryFile()
    f.write(vcf.encode()); f.flush()
    
    result = os.popen("cd /opt/vep/src/ensembl-vep ; ./vep -i %s -o STDOUT %s --assembly GRCh38  --offline --force_overwrite --dir /opt/vep/.vep/  " % (f.name,option) ).readlines()
    
    columns = ['Uploaded_variation','Location','Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation' ,'Extra']
    var_info =[]
    for x in result:
        if x[0] != '#':
            var =dict(zip(columns,x.split('\t') ) )
            extra       = var['Extra']
            # "a=b;c=d;e=f" => ['a=b','c=d','e=f']
            extra_items = extra.split(';')
            extra_dic = {}
            # ["a=b","c=d",..] => {"a":"b","c":"d",..}
            for i,item in enumerate(extra_items):
                key = item.split('=')[0]
                value = item.split('=')[1]
                # last_item  "a":"something\n" => "a":"something"
                if i ==len(extra_items)-1 and value[-1]=='n':
                    value = value[:-2]
                extra_dic[key]=value
            
            var['Extra'] = extra_dic
            var_info.append(var)    

    f.close()
    return var_info

def split_seq(seq):
    position = seq.split(' ')[0]
    variant  = seq.split(' ')[1]

    chrom    = position.split(':')[0]
    pos      = position.split(':')[1]

    ref      = variant.split('>')[0]
    alt      = variant.split('>')[1]
    return chrom , pos,ref,alt