# -*- coding: utf-8 -*-
__author__ = 'jiashun'
import gc

ch = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M']
for cc in ch:
    sam = open('/2_disk/xiaojs/sam/reads_align_'+str(cc)+'_sorted.sam', 'r')
    refs = open('/2_disk/longyk/human/redchr'+str(cc), 'r')
    pileup = open('/2_disk/xiaojs/pileup/pileup_'+str(cc), 'w')

    r0 = refs.readline()
    ref = refs.readline().strip()

    s0 = sam.readline()
    s0 = sam.readline()
    s0 = sam.readline()
    pile = ['' for i in xrange(len(ref)+1000)]
    qua = ['' for j in xrange(len(ref)+1000)]
    for s in sam:
        r = s.split('\t')
        if int(r[1]) == 0:
            if int(r[4]) == 42:
                lo = int(r[3])
                cigar = int(r[5].replace('M',''))
                quality = r[10]
                #print lo, cigar
                pile[lo-1] += '^K'
                for c in xrange(cigar):
                    pile[lo-1+c] += '.'
                    qua[lo-1+c] += quality[c]
                pile[lo-1+c] += '$'
                #print quality
            else:
                lo = int(r[3])
                seq = r[9]
                cigar = r[5].replace('M','').split('1X')
                quality = r[10]
                #print lo, cigar
                pile[lo-1] += '^I'
                for c in xrange(int(cigar[0])):
                    pile[lo-1+c] += '.'
                    qua[lo-1+c] += quality[c]
                pile[lo-1+int(cigar[0])] += seq[int(cigar[0])]
                qua[lo-1+int(cigar[0])] += quality[int(cigar[0])]
                for c in xrange(int(cigar[1])):
                    pile[lo+int(cigar[0])+c] += '.'
                    qua[lo+int(cigar[0])+c] += quality[int(cigar[0])+1+c]
                pile[lo+int(cigar[0])+c] += '$'
        elif int(r[1]) == 16:
            lo = int(r[3])
            quality = list(r[10].strip())
            quality.reverse()
            quality = ''.join(quality)
            if int(r[4]) == 42:
                cigar = int(r[5].replace('M',''))
                #print lo, cigar
                pile[lo-1] += '^K'
                for c in xrange(cigar):
                    pile[lo-1+c] += ','
                    qua[lo-1+c] += quality[c]
                pile[lo-1+c] += '$'
                #print quality
            else:
                seq = r[9].lower()
                cigar = r[5].replace('M','').split('1X')
                #print lo, cigar
                pile[lo-1] += '^I'
                for c in xrange(int(cigar[0])):
                    pile[lo-1+c] += ','
                    qua[lo-1+c] += quality[c]
                pile[lo-1+int(cigar[0])] += seq[int(cigar[0])]
                qua[lo-1+int(cigar[0])] += quality[int(cigar[0])]
                for c in xrange(int(cigar[1])):
                    pile[lo+int(cigar[0])+c] += ','
                    qua[lo+int(cigar[0])+c] += quality[int(cigar[0])+1+c]
                pile[lo+int(cigar[0])+c] += '$'

    #print pile

    for p in xrange(len(ref)):
        pileup.write('chromosome'+str(cc)+'\t'+str(p+1)+'\t'+ref[p]+'\t'+str(len(pile[p].replace('^K','').replace('^I','').replace('$','')))+'\t'+pile[p]+'\t'+qua[p]+'\n')

    sam.close()
    pileup.close()
    refs.close()
    gc.collect()