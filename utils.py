#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: utils.py
#          Desc: 辅助函数
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
    ref_range = {}
    for chrom in ref:
        avai_range = np.arange(0, len(ref[chrom][0]))
        n_base_indexes = re.finditer("[N]+", ref[chrom][0])
        for item in n_base_indexes:
            start = item.start(0)
            end = item.end(0)
            if not start - 20 <= 0:
                start = start - 20
            if not end + 20 >= len(ref[chrom][0]):
                end = end + 20
            item_range = np.arange(start, end)
            avai_range = np.delete(avai_range, item_range)

        ref_range[chrom] = avai_range
    return ref_range


def outputFa(ref, outfileName):
    outfile = open(outfileName, 'w')
    for key in ref:
        k = 0
        for line in ref[key]:
            i = 0
            k = k+1
            str_len = len(line)
            print str_len
            outfile.write('>'+key+'_'+str(k)+'\n')
            while i+50 <= str_len:
                outfile.write(line[i:i+50]+'\n')
                i = i+50
            if i < str_len:
                outfile.write(line[i:str_len]+'\n')
    outfile.close()


def read_pf_config(pf_file_name):
    pf_dict = {}
    with open(pf_file_name) as infile:
        for line in infile:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue
            listLine = line.split("\t")
            pf_dict[listLine[0]] = float(listLine[1])

    return pf_dict


def reference(ref_name):
    ref_dic = {}
    chr_name = ''
    tmp_str = ''
    for line in open(ref_name):
        newline = line.rstrip()
        if newline.startswith('>'):
            if '_' in newline:
                chr_split = newline.split('_')
            elif '-' in newline:
                chr_split = newline.split('-')
            else:
                chr_split = [newline]
            if chr_name != '':
                if chr_name not in ref_dic:
                    ref_dic[chr_name] = [tmp_str]
                else:
                    # 此处是含有paternal 和 maternal
                    ref_dic[chr_name].append(tmp_str)
            chr_name = chr_split[0].split('>')[1]
            tmp_str = ''
        else:
            tmp_str = tmp_str+newline
    if chr_name not in ref_dic:
        ref_dic[chr_name] = [tmp_str]
    else:
        ref_dic[chr_name].append(tmp_str)
    return ref_dic


def read_dbsnp(dbsnp, chome):
    snp_dic = {}
    if chome == 'ALL':
        for line in open(dbsnp):
            if not line.startswith('#'):
                newline = line.rstrip().split('\t')
                if not newline[0].startswith("chr"):
                    chromName = "chr{}".format(newline[0])
                else:
                    chromName = newline[0]

                if len(newline[4]) == 1 and len(newline[3]) == 1:
                    if chromName not in snp_dic:
                        snp_dic[chromName] = [
                            [newline[1], newline[3], newline[4]]]
                    else:
                        snp_dic[chromName].append(
                            [newline[1], newline[3], newline[4]])
    else:
        chrome_list = chome.split(',')
        print chrome_list
        for line in open(dbsnp):
            if not line.startswith('#'):
                newline = line.rstrip().split('\t')
                if not newline[0].startswith("chr"):
                    chromName = "chr{}".format(newline[0])
                else:
                    chromName = newline[0]

                if len(
                        newline[4]) == 1 and len(
                        newline[3]) == 1 and chromName in chrome_list:
                    if chromName not in snp_dic:
                        snp_dic[chromName] = [
                            [newline[1], newline[3], newline[4]]]
                    else:
                        snp_dic[chromName].append(
                            [newline[1], newline[3], newline[4]])
    for key in snp_dic.keys():
        print "key = {}".format(key)
        print "len(snp_dict[key]) = {}".format(len(snp_dic[key]))
    return snp_dic


def generate_normal(ref_dic, snp_dic, num, outfilename, hyp_rate=0.5):
    # 这个函数没有生成fasta文件
    outfile = open(outfilename, 'w')
    outfile.write('#chr\tpois\tref\talt\thap\n')
    snp_list = {}
    # hapl=[0,1]
    l1 = [1 for i in range(int(num*hyp_rate))]
    l2 = [0 for i in range(num-int(num*hyp_rate))]
    total_list = l1+l2
    random.shuffle(total_list)
    ref_list = []
    # 平均分配每一个染色体上？
    num = int(num/len(ref_dic.keys()))
    i = 0
    print "ref_dict keys():"
    print ref_dic.keys()
    for key in ref_dic.keys():
        print "ref_dic key ={}".format(key)
        print "ref len = {}".format(len(ref_dic[key]))
        print "snp len = {}".format(len(snp_dic[key]))
    print "snp_dict keys():"
    print snp_dic.keys()

    print "num = {}".format(num)

    for key in ref_dic:
        # 表示fasta 上是不是含有paternal 和maternal
        # 正常情况下hapl 是 [0]
        hapl = [k for k in range(len(ref_dic[key]))]
        tmp_str_list = ref_dic[key]
        str_list = []
        for tmp in tmp_str_list:
            str_list.append(list(tmp))
        if len(snp_dic[key]) >= num:
            random_list = random.sample(snp_dic[key], num)
            for line in random_list:
                if total_list[i] == 0:
                    # 随机获取一个位置
                    hap = random.sample(hapl, 1)[0]
                    # hap表示parternal or maternal ??
                    old = str_list[hap][int(line[0])-1]
                    str_list[hap][int(line[0])-1] = line[-1]
                    outfile.write(
                        key + '\t' + line[0] + '\t' + old + '\t' + line[-1] + '\t' +
                        str(hap + 1) + '\n')
                else:
                    for t in hapl:
                        old = str_list[t][int(line[0])-1]
                        str_list[t][int(line[0])-1] = line[-1]
                    outfile.write(
                        key + '\t' + line[0] + '\t' + old + '\t' + line[-1] +
                        '\thomozygous\n')
                i = i+1
                if key not in snp_list:
                    snp_list[key] = [int(line[0])]
                else:
                    snp_list[key].append(int(line[0]))
        else:
            for line in snp_dic[key]:
                if total_list[i] == 0:
                    hap = random.sample(hapl, 1)[0]
                    old = str_list[hap][int(line[0])-1]
                    str_list[hap][int(line[0])-1] = line[-1]
                    outfile.write(
                        key + '\t' + line[0] + '\t' + old + '\t' + line[-1] + '\t' +
                        str(hap + 1) + '\n')
                else:
                    for t in hapl:
                        old = str_list[t][int(line[0])-1]
                        str_list[t][int(line[0])-1] = line[-1]
                    outfile.write(
                        key + '\t' + line[0] + '\t' + old + '\t' + line[-1] +
                        '\thomozygous\n')
                i = i+1
                if key not in snp_list:
                    snp_list[key] = [int(line[0])]
                else:
                    snp_list[key].append(int(line[0]))
        tmp_str = []
        for tmp in str_list:
            tmp_str.append(''.join(tmp))
            # 将加入snp之后的放入ref_dic中，现在ref中的内容是另一条更改之后的，
            # 所以这里就没有paternal和maternal 之分了。
        ref_dic[key] = tmp_str
    outfile.close()
    return [ref_dic, snp_list]

