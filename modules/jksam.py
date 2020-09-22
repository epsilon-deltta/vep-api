#!/usr/bin/python

import sys, re
import jkgenome, jkbasic, jkbio


pattern_5p_softClip = re.compile('^([0-9]+)S')
pattern_3p_softClip = re.compile('([0-9]+)S$')
pattern_softClip = re.compile('([0-9]+)S')
pattern_hardClip = re.compile('([0-9]+)H')
pattern_splice = re.compile('(?=([A-Z]([0-9]+)M([0-9]+)N([0-9]+)M))')
#pattern_splice = re.compile('([0-9]+)N')
pattern_match = re.compile('([0-9]+)M')
pattern_cigar = re.compile('([0-9]+)([MIDNSHP])')
#pattern_anchor5 = re.compile('^([0-9]+)M[0-9]+I')
#pattern_anchor3 = re.compile('[0-9]+I([0-9]+)M$')

def processSAMFlag(flag):

    h = {}

    if flag & 0x0001:
        h['paired'] = True
    else:
        h['paired'] = False

    if flag & 0x0002:
        h['mapped_in_pair'] = True
    else:
        h['mapped_in_pair'] = False

    if flag & 0x0004:
        h['unmapped'] = True
    else:
        h['unmapped'] = False

    if flag & 0x0008:
        h['unmapped_mate'] = True
    else:
        h['unmapped_mate'] = False

    if flag & 0x0010:
        h['strand'] = '-'
    else:
        h['strand'] = '+'

    if flag & 0x0020:
        h['strand_mate'] = '-'
    else:
        h['strand_mate'] = '+'

    if flag & 0x0040:
        h['first_read_of_pair'] = True
    else:
        h['first_read_of_pair'] = False

    if flag & 0x0080:
        h['second_read_of_pair'] = True
    else:
        h['second_read_of_pair'] = False

    if flag & 0x0100:
        h['alignment_not_primary'] = True
    else:
        h['alignment_not_primary'] = False

    if flag & 0x0200:
        h['QC_failed'] = True
    else:
        h['QC_failed'] = False

    if flag & 0x0400:
        h['duplicate'] = True
    else:
        h['duplicate'] = False

    return h


class samLine:

    def __init__(self,line):

        tokL = line.rstrip().split('\t')

        self.readId = tokL[0]
        self.chrom = tokL[2]
        self.seq = tokL[9]
        self.pos = int(tokL[3])
        self.insertSize = int(tokL[8])
        self.cigar = tokL[5]
        self.flagH = processSAMFlag(int(tokL[1]))
        self.softclipping = sum([int(mo.group(1)) for mo in pattern_softClip.finditer(self.cigar)])
        self.hardclipping = sum([int(mo.group(1)) for mo in pattern_hardClip.finditer(self.cigar)])
        self.softclipping5p = sum([int(n) for n in pattern_5p_softClip.findall(self.cigar)])
        self.softclipping3p = sum([int(n) for n in pattern_3p_softClip.findall(self.cigar)])
        self.raw = line
        self.qual = tokL[10]

        for t in tokL[11:]:

            t = t.upper()

            if t[:4] == 'NH:I':
                self.NH = int(t[5:])

            if t[:4] == 'NM:I':
                self.NM = int(t[5:])

            if t[:4] == 'XM:I':
                self.XM = int(t[5:])

            if t[:4] == 'MD:Z':
                self.MD = t[5:]

            # HISAT
            if t[:4] == 'XS:A':
                self.trans_strand = t[5:]

            # HISAT-specific
            if t[:4] == 'YT:Z':
                self.YT = t[5:]
#
#
#            # STAR-specific
#            if t[:4] == 'jM:B':
#                self.jM = t[7:].split(',')

        if 'NH:i:' not in line:
            self.NH = 1
            self.raw = '%s\tNH:i:1\n' % (self.raw[:-1],)

    def mateId(self):
        return min(self.seq,jkbasic.rc(self.seq))

#    def fastq(self):
#
#        if self.flagH['strand'] == '-':
#            seq_ori = jkbio.rc(self.seq)
#            qual_ori = jkbio.rev(self.qual)
#        else:
#            seq_ori = self.seq
#            qual_ori = self.qual
#        
#        return '@%s\n%s\n+\n%s\n' % (self.readId,seq_ori,qual_ori)

#    def containsCorrectSplice(self,minIntron=20,maxIntron=500000,minAnchor=3): 
#
#        if not 'N' in self.cigar:
#            return False
#
#        flag = False
#
#        mo_groupL = pattern_splice.findall('X'+self.cigar)
#
##        if self.cigar.count('N') != len(mo_groupL):
##            print 'junction containing indel:', self.cigar, self.cigar.count('N'), mo_groupL, len(mo_groupL)
#
#        for mo_g in mo_groupL:
#
#            left_anchor = int(mo_g[1])
#            intron = int(mo_g[2])
#            right_anchor = int(mo_g[3])
#
#            if left_anchor >= minAnchor and right_anchor >= minAnchor and minIntron <= intron <= maxIntron:
#                flag = True
#                break
#
#        return flag

    def containsCorrectSplice(self,minIntron=20,maxIntron=500000,minAnchor=3): 

        if not 'N' in self.cigar:
            return False

        flag = False

        mo_groupL = pattern_cigar.findall(self.cigar)

        for i in range(len(mo_groupL)):

            mo_g = mo_groupL[i]

            bases = int(mo_g[0])
            alt = mo_g[1]

            if alt == 'N':

                left_bases = int(mo_groupL[i-1][0])
                left_alt = mo_groupL[i-1][1]
                right_bases = int(mo_groupL[i+1][0])
                right_alt = mo_groupL[i+1][1]

                if minIntron <= bases <= maxIntron and left_alt == 'M' and left_bases >= minAnchor and right_alt == 'M' and right_bases >= minAnchor:
                    flag = True
                    break

        return flag


#    def spliceMax(self): # zero if not cover a splice junction
#
#        spliceMax = 0
#        
#        for mo in pattern_splice.finditer(self.cigar):
#
#            spliceSize = int(mo.group(1))
#            if spliceSize > spliceMax:
#                spliceMax = spliceSize
#
#        return spliceMax

#    def minAnchorSize(self): # zero if not cover a splice junction
#
#        anchorMin = 1000
#        
#        mo = pattern_anchor5.match(self.cigar)
#
#        if mo and int(mo.group(1)) < anchorMin:
#            anchorMin = int(mo.group(1))
#    
#        mo = pattern_anchor3.match(self.cigar)
#
#        if mo and int(mo.group(1)) < anchorMin:
#            anchorMin = int(mo.group(1))
# 
#        if anchorMin == 1000:
#            return 0
#        else:
#            return anchorMin

    def getSplicingCoordinates(self,minIntron=20,maxIntron=500000,minAnchor=3,minReadLen=30): # first and last base (1-base) of intron

        spliceL = []
        curChrPos = self.pos

        if len(self.seq) < minReadLen or self.softclipping5p > 0:
            return []

        if hasattr(self,'trans_strand') and self.trans_strand in ('+','-'):
            strand = self.trans_strand
#            if self.trans_strand == self.flagH['strand']:
#                strand = '+'
#            else:
#                strand = '-'
        else:
            strand = '?'

        mo_groupL = pattern_cigar.findall(self.cigar)

        for i in range(len(mo_groupL)):

            mo_g = mo_groupL[i]

            bases = int(mo_g[0])
            alteration = mo_g[1]

            if alteration == 'P':
                print(self.cigar)
                raise Exception

            if alteration in ('H','S','I'):
                continue

            if alteration == 'N':

                left_bases = int(mo_groupL[i-1][0])
                left_alt = mo_groupL[i-1][1]
                right_bases = int(mo_groupL[i+1][0])
                right_alt = mo_groupL[i+1][1]

                if left_alt == 'M' and right_alt == 'M' and not (minIntron <= bases <= maxIntron):
                    return []

                if left_alt == 'M' and left_bases >= minAnchor and right_alt == 'M' and right_bases >= minAnchor:
                    spliceL.append((self.chrom,curChrPos,curChrPos+bases-1,strand))

            curChrPos += bases

        return spliceL

    def matchL(self):
        return sum([int(mo.group(1)) for mo in pattern_match.finditer(self.cigar)])

    def raw_multimapStripInfo(self):
        
        tL = self.raw.split('\t')

        tL[5] = '*'
        tL[11] = 'NH:i:%s' % (self.NH,)

        return '\t'.join(tL[:12])


class samFile(file):

    def __init__(self, fileName, unit='block'): # unit is block (lines with same consecutive IDs) or line (1-line)

        f = super(samFile,self).__init__(fileName)

        self.unit = unit

        headerL = []

        while 1:

            try:
                self.firstLine = file.next(self)
            except StopIteration:
                self.firstLine = None
                break

            if self.firstLine[0] =='@':
                headerL.append(self.firstLine)
            else:
                break

        self.header = ''.join(headerL)

        return f
    
    def next(self):

        if self.unit=='block':

            if self.firstLine == None:
                raise StopIteration

            samL = [samLine(self.firstLine)]
            self.firstLineSeqID = self.firstLine.split('\t')[0]
                            
            while 1:

                try:
                    nextLine = file.next(self)
                except StopIteration:
                    self.firstLine = None
                    break
                
                if nextLine.split('\t')[0] != self.firstLineSeqID: 
                    # this is the last record or the last line of this record
                    self.firstLine = nextLine
                    break
                else:
                    samL.append(samLine(nextLine))

            return samL

        elif self.unit=='line':

            if self.firstLine == None:
                raise StopIteration

            sam = samLine(self.firstLine)

            try:
                self.firstLine = file.next(self)
            except StopIteration:
                self.firstLine = None
            
            return sam

        else:

            raise Exception
