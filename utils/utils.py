#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: utils.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-07-26 13:46:41
#       History:
# =============================================================================
'''

import numpy as np
import re
import random


def compute_range(ref):
    refRange = {}
    for chrom in ref:
        avaiRange = np.arange(0, len(ref[chrom]['P'][0]))
        nBaseIndexes = re.finditer("[N]+", ref[chrom]['P'][0])
        for item in nBaseIndexes:
            start = item.start(0)
            end = item.end(0)
            if not start - 20 <= 0:
                start = start - 20
            if not end + 20 >= len(ref[chrom]['P'][0]):
                end = end + 20
            itemRange = np.arange(start, end)
            avaiRange = np.delete(avaiRange, itemRange)

        refRange[chrom] = avaiRange
    return refRange


def outputFa(ref, outfileName):
    outfile = open(outfileName, 'w')

    for chrom in ref.keys():
        for haplType in ref[chrom].keys():
            for haplIdx in range(len(ref[chrom][haplType])):
                strID = '>{0}_{1}_{2}\n'.format(chrom, haplType, haplIdx)
                outfile.write(strID)

                strDNA = ref[chrom][haplType][haplIdx]
                strDNALen = len(strDNA)
                i = 0

                while i+50 <= strDNALen:
                    outfile.write(strDNA[i:i+50]+'\n')
                    i = i+50

                if i < strDNALen:
                    outfile.write(strDNA[i:strDNALen]+'\n')

    outfile.close()


def read_pf_config(pf_file_name):
    pfDict = {}
    with open(pf_file_name) as infile:
        for line in infile:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue
            listLine = line.split("\t")
            pfDict[listLine[0]] = float(listLine[1])

    return pfDict


def reference(ref_name):
    refDic = {}
    chrName = ''
    tmpStr = ''

    for line in open(ref_name):
        newline = line.strip()

        if newline.startswith('>'):
            if '-' in newline:
                chrSplit = newline.split('-')
            else:
                chrSplit = [newline]

            if chrName != '':
                if chrName not in refDic:
                    refDic[chrName] = {"P": [tmpStr]}
                else:
                    refDic[chrName]["M"] = [tmpStr]

            chrName = chrSplit[0].split('>')[1]
            tmpStr = ''
        else:
            tmpStr = tmpStr+newline

    if chrName not in refDic:
        refDic[chrName] = {"P": [tmpStr]}
    else:
        refDic[chrName]["M"] = [tmpStr]
    return refDic


def read_dbsnp(dbsnp, chome):
    snpDic = {}
    if chome == 'ALL':
        for line in open(dbsnp):
            if not line.startswith('#'):
                newline = line.rstrip().split('\t')
                if not newline[0].startswith("chr"):
                    chromName = "chr{}".format(newline[0])
                else:
                    chromName = newline[0]

                if len(newline[4]) == 1 and len(newline[3]) == 1:
                    if chromName not in snpDic:
                        snpDic[chromName] = [
                            [newline[1], newline[3], newline[4]]]
                    else:
                        snpDic[chromName].append(
                            [newline[1], newline[3], newline[4]])
    else:
        chromeList = chome.split(',')
        for line in open(dbsnp):
            if not line.startswith('#'):
                newline = line.rstrip().split('\t')
                if not newline[0].startswith("chr"):
                    chromName = "chr{}".format(newline[0])
                else:
                    chromName = newline[0]

                if len(
                        newline[4]) == 1 and len(
                        newline[3]) == 1 and chromName in chromeList:
                    if chromName not in snpDic:
                        snpDic[chromName] = [
                            [newline[1], newline[3], newline[4]]]
                    else:
                        snpDic[chromName].append(
                            [newline[1], newline[3], newline[4]])
    return snpDic


def generate_normal(refDic, snpDic, num, outfilename, hyp_rate=0.5):
    # 这个函数没有生成fasta文件
    outfile = open(outfilename, 'w')
    outfile.write('#chr\tpois\tref\talt\thap\n')
    snpList = {}
    # hapl=[0,1]
    l1 = [1 for i in range(int(num*hyp_rate))]
    l2 = [0 for i in range(num-int(num*hyp_rate))]
    totalList = l1+l2
    random.shuffle(totalList)
    # 平均分配每一个染色体上？
    num = int(num/len(refDic.keys()))
    i = 0

    for key in refDic:
        # 表示fasta 上是不是含有paternal 和maternal
        # 正常情况下hapl 是 [0]
        # hapl = [k for k in range(len(refDic[key]))]
        hapl = refDic[key].keys()

        tmpStrList = [refDic[key][hapl[0]][0], refDic[key][hapl[1]][0]]

        strList = []
        for tmp in tmpStrList:
            strList.append(list(tmp))
        if len(snpDic[key]) >= num:
            randomList = random.sample(snpDic[key], num)
            for line in randomList:
                if totalList[i] == 0:
                    # 随机获取一个位置
                    hap = random.sample(range(len(hapl)), 1)[0]
                    # hap表示parternal or maternal ??
                    old = strList[hap][int(line[0])-1]
                    strList[hap][int(line[0])-1] = line[-1]
                    outfile.write(
                        key + '\t' + line[0] + '\t' + old + '\t' + line[-1] + '\t' +
                        str(hap + 1) + '\n')
                else:
                    for t in range(len(hapl)):
                        old = strList[t][int(line[0])-1]
                        strList[t][int(line[0])-1] = line[-1]
                    outfile.write(
                        key + '\t' + line[0] + '\t' + old + '\t' + line[-1] +
                        '\thomozygous\n')
                i = i+1
                if key not in snpList:
                    snpList[key] = [int(line[0])]
                else:
                    snpList[key].append(int(line[0]))
        else:
            for line in snpDic[key]:
                if totalList[i] == 0:
                    hap = random.sample(range(len(hapl)), 1)[0]
                    old = strList[hap][int(line[0])-1]
                    strList[hap][int(line[0])-1] = line[-1]
                    outfile.write(
                        key + '\t' + line[0] + '\t' + old + '\t' + line[-1] + '\t' +
                        str(hap + 1) + '\n')
                else:
                    for t in range(len(hapl)):
                        old = strList[t][int(line[0])-1]
                        strList[t][int(line[0])-1] = line[-1]
                    outfile.write(
                        key + '\t' + line[0] + '\t' + old + '\t' + line[-1] +
                        '\thomozygous\n')
                i = i+1
                if key not in snpList:
                    snpList[key] = [int(line[0])]
                else:
                    snpList[key].append(int(line[0]))
        tmpStr = []
        for i in range(len(strList)):
            tmpStr = ''.join(strList[i])
            # 将加入snp之后的放入refDic中，现在ref中的内容是另一条更改之后的，
            # 所以这里就没有paternal和maternal 之分了。
            refDic[key][hapl[i]] = [tmpStr]
    outfile.close()
    return [refDic, snpList]


def rand_DNA(length=10):
    return ''.join(random.choice(['A', 'T', 'G', 'C']) for _ in range(length))


def SNP(str1):
    l = ['A', 'T', 'C', 'G', 'N']
    l.remove(str1)
    while True:
        tmp = random.sample(l, 1)[0]
        if tmp != 'N':
            break
        else:
            continue
    return tmp
