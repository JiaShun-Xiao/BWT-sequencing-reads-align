# -*- coding: utf-8 -*-
__author__ = 'jiashun'
import re
import gc
depth = 6
threshold = 0.4
asc = {'z': 122, 'q': 113, "}": 125, 'p': 112, 'h': 104, 'D': 68, 'X': 88, 'Y': 89, 'o': 111, 'A': 65, 'B': 66,
       'C': 67, 'R': 82, 'S': 83, 'P': 80, 'Q': 81, 'V': 86, 'W': 87, 'T': 84, 'U': 85, '{': 123, 'e': 101, '_': 95,
       'f': 102, 'i': 105, 'k': 107, 't': 116, 'a': 97, '`': 96, 's': 115, 'x': 120, 'd': 100, 'w': 119, 'M': 77,
       'L': 76, 'K': 75, 'J': 74, 'I': 73, 'H': 72, 'G': 71, 'F': 70, '[': 91, 'Z': 90, ']': 93,
       '^': 94, 'O': 79, 'N': 78, 'v': 118, 'c': 99, 'y': 121, 'b': 98, 'l': 108, 'u': 117, 'm': 109, 'j': 106,
       'r': 114, '|': 124, 'g': 103, 'n': 110, 'E': 69, '\\': 92, '%': 37, '2': 50, '<': 60, '=': 61, '>': 62, '?': 63,
       '@': 64, '5': 53, '4': 52, '0': 48, '1': 49, '.': 46, '/': 47, ',': 44, '-': 45, '*': 42, '+': 43, '(': 40,
       ')': 41, '\'': 39, '&': 38, ';': 59, ':': 58, '!': 33, '6': 54, '9': 57, '8': 56, '3': 51, '$': 36, '#': 35,
        '"': 34, '7': 55}


def compare(a,t,g,c):
    if max(a,t,g,c) == a:
        return 'A'
    elif max(a,t,g,c) == t:
        return 'T'
    elif max(a,t,g,c) == g:
        return 'G'
    elif max(a,t,g,c) == c:
        return 'C'


ch = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M']
for cc in ch:
    pileup = open('/2_disk/xiaojs/pileup/pileup_'+str(cc), 'r')
    vcf = open('/2_disk/xiaojs/vcf/vcf_'+str(cc), 'w')
    vcf.write('Min coverage:'+'\t'+str(depth)+'\n')
    vcf.write('Min var freq:'+'\t'+str(threshold)+'\n')
    vcf.write('Chrom'+str(cc)+'\t'+'Position'+'\t'+'Ref'+'\t'+'Var'+'\t'+'Qual'+'\t'+'Cons:Cov:Reads1:Reads2:Freq:P-value'+'\t'+'StrandFilter:R1+:R1-:R2+:R2-:pval'+'\n')

    for line in pileup:
        p = line.strip().split('\t')
        #print p
        if len(p) > 4:
            seq = p[4].replace('^K','').replace('^I','').replace('$','')
            if re.match(r'(.*)[ATGCatgc]+(.*)', seq) and int(p[3]) >= depth:
                #print seq
                A = 0
                a = 0
                T = 0
                t = 0
                G = 0
                g = 0
                C = 0
                c = 0
                dot = 0
                for li in range(len(seq)):
                    if seq[li] == 'A':
                        A += 1
                    elif seq[li] == 'a':
                        a += 1
                    elif seq[li] == 'T':
                        T += 1
                    elif seq[li] == 't':
                        t += 1
                    elif seq[li] == 'G':
                        G += 1
                    elif seq[li] == 'g':
                        g += 1
                    elif seq[li] == 'C':
                        C += 1
                    elif seq[li] == 'c':
                        c += 1
                    else:
                        dot += 1
                frequence = (A+a+T+t+G+g+C+c)/float(dot+A+a+T+t+G+g+C+c)
                if dot == 0 or frequence >= threshold:
                    #print seq
                    qual = p[5]
                    #print qual
                    quals = 0
                    for q in qual:
                        quals += int(asc[q])
                    quals /= len(qual)
                    vcf.write(p[0]+'\t'+p[1]+'\t'+p[2]+'\t'+compare(A+a,T+t,G+g,C+c)+'\t'+str(quals)+'\t'+
                              compare(A+a,T+t,G+g,C+c)+':'+str(len(qual))+':0:'+str(A+a+T+t+G+g+C+c)+':'+
                              str(frequence)+':0'+'\t'+'Pass:0:0:'+str(A+T+C+G)+':'+str(a+g+t+c)+':1E0'+'\n')

    pileup.close()
    vcf.close()
    gc.collect()