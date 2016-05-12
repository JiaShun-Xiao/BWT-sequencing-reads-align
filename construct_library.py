# -*- coding: utf-8 -*-
__author__ = 'jiashun'
import gc

ch = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M']


def rotations(t):
    global tt
    rt = []
    for i in xrange(len(t)):
        rt.append(tt[i:i+co])
    return rt


def flat(seq):
    res = []
    for i in seq:
        if isinstance(i, list):
            res.extend(flat(i))
        else:
            res.append(i)
    return res

co = 50
for c in ch:
    print c
    op = 1
    #file = open('hba.txt','r')
    #file = open('lambda_virus.fa','r')
    #file = open('/2_disk/xiaojs/python/blast/chrom/chromosome14.txt','r')
    file = open('/2_disk/longyk/human/redchr'+str(c), 'r')
    first_colum = open('/2_disk/xiaojs/human/first_colum_'+str(c), 'w')
    bwt = open('/2_disk/xiaojs/human/bwt_'+str(c), 'w')
    suffix = open('/2_disk/xiaojs/human/suffix_'+str(c), 'w')
    tally = open('/2_disk/xiaojs/human/tally_'+str(c), 'w')

    t0 = file.readline()
    t = file.readline().strip()
    A = 0
    T = 0
    G = 0
    C = 0
    for li in range(len(t)):
        if t[li] == 'A':
            A += 1
        elif t[li] == 'T':
            T += 1
        elif t[li] == 'G':
            G += 1
        elif t[li] == 'C':
            C += 1
    first_colum.write(str(A)+"\n"+str(C)+"\n"+str(G)+"\n"+str(T))

    t = t+'$'
    tt = t*2
    last = tt[len(t)-1:len(tt)-1]
    rt = rotations(t)
    #find first 20 sequence and their corresponding location after rotation
    rt_s = {}
    #sort first 20 sequence after rotation
    rt_ss = []

    ops = False
    for i in xrange(0, len(rt)):
        if rt_s.has_key(rt[i]):
            #print 'no enough column!'
            #print rt[i][0:co]
            ops = True
            if not isinstance(rt_s[rt[i]], list):
                rt_s[rt[i]] = [rt_s[rt[i]]]
            rt_s[rt[i]].append(i)
        else:
            rt_s[rt[i]] = i
            rt_ss.append(rt[i])
    #print rt_s

    srt = sorted(rt_ss)
    #print srt

    while ops:
        op += 2
        print 'op', op
        ops = False
        srt = flat(srt)
        for ind, i in enumerate(srt):
            if isinstance(rt_s[i], list):
                #print rt_s[i]
                #print i,
                rt_ss2 = []
                for j in rt_s[i]:
                    if rt_s.has_key(tt[j:j+op*co]):
                        #print 'no enough column!'
                        #print tt[j:j+op*co]
                        #print rt_s[tt[j:j+op*co]]
                        if not isinstance(rt_s[tt[j:j+op*co]], list):
                            rt_s[tt[j:j+op*co]] = [rt_s[tt[j:j+op*co]]]
                        rt_s[tt[j:j+op*co]].append(j)
                        ops = True
                    else:
                        rt_s[tt[j:j+op*co]] = j
                        rt_ss2.append(tt[j:j+op*co])
                temp = sorted(rt_ss2)
                #print 'index', ind
                srt[ind] = temp

        #print rt_s
        #print srt
        #print ops

    #print srt
    srt = flat(srt)
    #print srt
    last_colum = ''
    for i in srt:
        last_colum += last[rt_s[i]]
        suffix.write(str(rt_s[i])+" ")
    bwt.write(last_colum)
    tallys = [ [0] for i in xrange(4)]
    for i in range(len(last_colum)):
        if last_colum[i] == 'A':
            tallys[0].append(1+tallys[0][-1])
            tallys[1].append(0+tallys[1][-1])
            tallys[2].append(0+tallys[2][-1])
            tallys[3].append(0+tallys[3][-1])
        elif last_colum[i] == 'C':
            tallys[0].append(0+tallys[0][-1])
            tallys[1].append(1+tallys[1][-1])
            tallys[2].append(0+tallys[2][-1])
            tallys[3].append(0+tallys[3][-1])
        elif last_colum[i] == 'G':
            tallys[0].append(0+tallys[0][-1])
            tallys[1].append(0+tallys[1][-1])
            tallys[2].append(1+tallys[2][-1])
            tallys[3].append(0+tallys[3][-1])
        elif last_colum[i] == 'T':
            tallys[0].append(0+tallys[0][-1])
            tallys[1].append(0+tallys[1][-1])
            tallys[2].append(0+tallys[2][-1])
            tallys[3].append(1+tallys[3][-1])
        elif last_colum[i] == '$':
            tallys[0].append(0+tallys[0][-1])
            tallys[1].append(0+tallys[1][-1])
            tallys[2].append(0+tallys[2][-1])
            tallys[3].append(0+tallys[3][-1])
    for i in xrange(4):
        del tallys[i][0]
    for i in xrange(4):
        for j in xrange(len(tallys[i])):
            tally.write(str(tallys[i][j])+" ")
        tally.write("\n")
    file.close()
    bwt.close()
    first_colum.close()
    tally.close()
    suffix.close()
    gc.collect()
