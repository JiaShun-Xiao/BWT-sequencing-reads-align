# -*- coding: utf-8 -*-
__author__ = 'jiashun'
import gc

ch = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','M']
for cc in ch:
    vcf = open('/2_disk/xiaojs/vcf/vcf_'+str(cc), 'r')
    vcf_new = open('/2_disk/xiaojs/vcf_new/vcf_new_'+str(cc), 'w')
    lo = open('/2_disk/longyk/human/locachr'+str(cc), 'r')
    los = []
    for l in lo:
        los.append(int(l))
    v0 = vcf.readline()
    v1 = vcf.readline()
    v2 = vcf.readline()
    vcf_new.write(v0)
    vcf_new.write(v1)
    vcf_new.write(v2)
    for line in vcf:
        v = line.strip().split('\t')
        loc = int(v[1])
        v[1] = str(los[loc/50] + loc%50)
        v = '\t'.join(v)
        vcf_new.write(v+'\n')
