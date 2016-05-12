# -*- coding: utf-8 -*-
__author__ = 'jiashun'
import re
from operator import itemgetter
import gc



def align(query):
    if query == '':
        return -1
    que = list(query)
    que.reverse()
    pre = [0, 0]
    sur = [0, 0]
    le = len(que)
    for q in xrange(le):
        if que[q] == 'N':
            ques = que[q+1:le]
            ques.reverse()
            ques = ''.join(ques)
            return align(ques)
        elif que[q] == 'A':
            if q == 0:
                pre = [1, A]
            else:
                if bwt[pre[0]] == 'A':
                    sur = [tallys[0][pre[0]],tallys[0][pre[1]]]
                    pre = [tallys[0][pre[0]],tallys[0][pre[1]]]
                else:
                    sur = [tallys[0][pre[0]]+1,tallys[0][pre[1]]]
                    pre = [tallys[0][pre[0]]+1,tallys[0][pre[1]]]
        elif que[q] == 'C':
            if q == 0:
                pre = [A+1,A+C]
            else:
                if bwt[pre[0]] == 'C':
                    sur = [tallys[1][pre[0]],tallys[1][pre[1]]]
                    pre = [A+tallys[1][pre[0]],A+tallys[1][pre[1]]]
                else:
                    sur = [tallys[1][pre[0]]+1,tallys[1][pre[1]]]
                    pre = [A+tallys[1][pre[0]]+1,A+tallys[1][pre[1]]]
        elif que[q] == 'G':
            if q == 0:
                pre = [A+C+1,A+C+G]
            else:
                if bwt[pre[0]] == 'G':
                    sur = [tallys[2][pre[0]],tallys[2][pre[1]]]
                    pre = [A+C+tallys[2][pre[0]],A+C+tallys[2][pre[1]]]
                else:
                    sur = [tallys[2][pre[0]]+1,tallys[2][pre[1]]]
                    pre = [A+C+tallys[2][pre[0]]+1,A+C+tallys[2][pre[1]]]
        elif que[q] == 'T':
            if q == 0:
                pre = [A+C+G+1,A+C+G+T]
            else:
                if bwt[pre[0]] == 'T':
                    sur = [tallys[3][pre[0]],tallys[3][pre[1]]]
                    pre = [A+C+G+tallys[3][pre[0]],A+C+G+tallys[3][pre[1]]]
                else:
                    sur = [tallys[3][pre[0]]+1,tallys[3][pre[1]]]
                    pre = [A+C+G+tallys[3][pre[0]]+1,A+C+G+tallys[3][pre[1]]]
        if sur[0]>sur[1]:
            return -1

    for i in range(sur[0], sur[1]+1):
        if que[le-1] == 'A':
            return suffix[i]+1
        elif que[le-1] == 'C':
            return suffix[A+i]+1
        elif que[le-1] == 'G':
            return suffix[A+C+i]+1
        elif que[le-1] == 'T':
            return suffix[A+C+G+i]+1


def complement(query):
    que_ = list(query)
    que_.reverse()
    com = []
    for q in que_:
        if q == 'A':
            com.append('T')
        elif q == 'C':
            com.append('G')
        elif q == 'G':
            com.append('C')
        elif q == 'T':
            com.append('A')
        else:
            com.append(q)
    return ''.join(com)


def alignment(qu, re):
    mm = 0
    temp = 0
    for i in xrange(len(re)):
        if qu[i] == 'N':
            continue
        else:
            if qu[i] != re[i]:
                mm += 1
                temp = i
    if mm == 1:
        return temp
    else:
        return -1


def mismatch(seq):
    seq_1 = seq[0:len(seq)/2]
    seq_2 = seq[len(seq)/2:len(seq)]
    le = len(seq_1)
    #print seq_1,seq_2
    global lo, flag, mapq, cigar
    seq_1 = list(seq_1)
    while seq_1[-1] == 'N':
        #print seq_1[-1]
        seq_1.pop()
    seq_1 = ''.join(seq_1)
    lo = align(seq_1)
    if lo > 0:
        #print lo
        m = alignment(seq, ref[lo-1:lo+len(seq)-1])
        #print seq,ref[lo-1:lo+len(seq)-1]
        #print m
        if m > 0:
            flag = 0
            mapq = 40
            cigar = str(m)+'M'+'1X'+str(len(seq)-m-1)+'M'
        else:
            flag = 4
            mapq = 0
            cigar = '*'
            lo = -1
    else:
        lo = align(seq_2)
        #print lo
        if lo > 0:
            m = alignment(seq, ref[lo-len(seq)/2-1:lo+len(seq)/2-1])
            #print seq
            #print ref[lo-len(seq)/2-1:lo+len(seq)/2-1]
            #print m
            if m > 0:
                lo -= le
                flag = 0
                mapq = 40
                cigar = str(m)+'M'+'1X'+str(len(seq)-m-1)+'M'
            else:
                flag = 4
                mapq = 0
                cigar = '*'
                lo = -1


ch = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M']
for c in ch:
    first_colum = open('/2_disk/xiaojs/human/first_colum_'+str(c), 'r')
    bwts = open('/2_disk/xiaojs/human/bwt_'+str(c), 'r')
    suffixs = open('/2_disk/xiaojs/human/suffix_'+str(c), 'r')
    tally = open('/2_disk/xiaojs/human/tally_'+str(c), 'r')
    sam = open('/2_disk/xiaojs/sam/reads_align_'+str(c)+'_sorted.sam', 'w')
    reads = open('/1_disk/public_resources/jiankuihe-exome-1.fq', 'r')
    reads_r = open('/1_disk/public_resources/jiankuihe-exome-1.fq', 'r')
    #reads = open('test_xiao.fq', 'r')
    #reads_r = open('test_xiao.fq', 'r')
    #refs = open('/2_disk/xiaojs/python/blast/chrom/chromosome14.txt', 'r')
    refs = open('/2_disk/longyk/human/redchr'+str(c), 'r')
    result = open('result.txt', 'a')

    r0 = refs.readline()
    ref = refs.readline().strip()

    A = int(first_colum.readline().strip())
    C = int(first_colum.readline().strip())
    G = int(first_colum.readline().strip())
    T = int(first_colum.readline().strip())

    bwt = bwts.read().strip()

    sam.write('@HD'+'\t'+'VN:1.0'+'\t'+'SO:sorted'+'\n')
    sam.write('@SQ'+'\t'+'chromosome'+str(c)+'\t'+'LN:'+str(len(bwt)-1)+'\n')
    sam.write('@PG'+'\t'+'ID:xiaojs'+'\t'+'PN:xiaojs'+'\t'+'VN:1.0'+'\n')

    suffix = suffixs.read().strip()
    suffix = suffix.split()
    suffix = map(int, suffix)

    tallys = [[0] for i in xrange(4)]
    temp = tally.readline().strip()
    temp = temp.split()
    temp = [int(i) for i in temp]
    tallys[0] = temp

    temp = tally.readline().strip()
    temp = temp.split()
    temp = [int(i) for i in temp]
    tallys[1] = temp

    temp = tally.readline().strip()
    temp = temp.split()
    temp = [int(i) for i in temp]
    tallys[2] = temp

    temp = tally.readline().strip()
    temp = temp.split()
    temp = [int(i) for i in temp]
    tallys[3] = temp
    r = 0
    ri = 0
    sort = {}
    flag = 0
    mapq = 0
    cigar = ''
    for line in reads:
        if re.match(r"@(.*)", line):
            reads_r.readline().strip()
            r += 1
            #print 'r  ', r
            seq = reads_r.readline().strip()
            reads_r.readline()
            qual = reads_r.readline().strip()
            if seq == '':
                continue
            seq = list(seq)
            qual = list(qual)
            while seq[0] == 'N':
                seq.pop(0)
                qual.pop(0)
            seq.reverse()
            qual.reverse()
            while seq[0] == 'N':
                seq.pop(0)
                qual.pop(0)
            seq.reverse()
            qual.reverse()
            seq = ''.join(seq)
            qual = ''.join(qual)
            #print r, seq, qual
            lo = align(seq)
            if lo > 0:
                flag = 0
                mapq = 42
                cigar = str(len(seq))+'M'
            else:
                seq_r = complement(seq)
                lo = align(seq_r)
                if lo > 0:
                    flag = 16
                    mapq = 42
                    cigar = str(len(seq))+'M'
                else:
                    mismatch(seq)
                    if lo < 0:
                        #print seq_r
                        mismatch(seq_r)
                        if flag == 0:
                            flag = 16
            if lo < 0:
                sort['r'+str(r)+'\t'+str(flag)+'\t'+'*'+'\t'+'0'+'\t'+str(mapq)+'\t'+cigar+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq+'\t'+qual+'\n'] = -1
            else:
                ri += 1
                #print 'ri ', ri
                #print 'r'+str(r)+'\t'+str(flag)+'\t'+'chromosome1'+'\t'+str(lo)+'\t'+str(mapq)+'\t'+cigar+'\n'
                sort['r'+str(r)+'\t'+str(flag)+'\t'+'chromosome'+str(c)+'\t'+str(lo)+'\t'+str(mapq)+'\t'+cigar+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq+'\t'+qual+'\n'] = lo
            #print seq
            #print complement(seq)
    #print sort
    sort = sorted(sort.iteritems(), key=itemgetter(1), reverse=False)
    #print sort
    for s in sort:
        if s[1] > 0:
            sam.write(s[0])
    result.write('chromosome'+str(c)+'\n')
    result.write('total '+str(r)+' reads'+'\n')
    result.write(str(ri)+' reads were successful aligned'+'\n')


    bwts.close()
    first_colum.close()
    tally.close()
    suffixs.close()
    reads.close()
    reads_r.close()
    result.close()
    sam.close()
    refs.close()
    gc.collect()