#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: breakpoint.py
#          Desc: breakpoints
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-07-28 10:10:00
#       History:
# =============================================================================
'''
from collections import Counter
import random
import copy


class Break:

    def __init__(self):
        self.chrom = ""
        self.hapl_type = -1
        self.hapl_index = -1
        self.position = -1
        self.name = ""
        self.insertStr = None

        # deletion的paired position
        # l
        self.paired_position = -1


class BreakPoints():
    """Generate BreakPoints for SCNA
   注意：
   1）生成breakpoint位置时候不应与已有变异重合，即应该在available空间之中
   2）生成breakpoint位置与available等现场，可以保存
   3）可以额外添加变异，恢复现场，生成breakpoint
    """

    def __init__(self):
        """Initialize from the variant dict."""

        # 原始ref，用来生成变异序列
        # 返回值不能用ref做为应用，确保ref引用始终是原始，即不能向ref赋值
        # self._ref = ref
        # available 空间，这个空间应该是父节点去除了变异位置之后的空间
        # available候选空间在breakpoint和变异位置生成时要兼容，以防止变异重叠
        # 这样插入位置前后1个bp都要去除
        # 这个变量随时改变，要向该变量赋值, 作为传递变量
        # self._ref_avail = ref_avail
        # 当前节点的变异内容，只需要额外添加该变异
        # dic value: 位置、变异长度、拷贝数、基因型、名称

        # 额外添加变异， BreakPoints 从父节点继承
        # self._variant_dict = variant_dict

        # 现场信息
        # 生成变异区间，用来生成变异位置时候，判断并且确保不出现插入重合现象
        self.breaks_list = []

    def generateBPs(self, variant_dict, avail_position, ref):
        self._generateDupBPs(variant_dict, avail_position, ref)
        self._generateDelBPs(variant_dict)

    def generateFa(self, ref):
        newFa = copy.deepcopy(ref)
        for chrom in ref.keys():
            if not self._bpsHasChrom(chrom):
                continue
            else:
                for hap_i in range(len(ref[chrom])):
                    if not self._bpsHasChromHap(chrom, hap_i):
                        continue
                    else:
                        newFa[chrom][hap_i] = self._generateHapStr(chrom, hap_i,
                                                                   ref)
        return newFa

    def _bpsHasChrom(self, chrom):
        bks = filter(lambda item: item.chrom == chrom, self.breaks_list)
        if len(bks) > 0:
            return True
        else:
            return False

    def _bpsHasChromHap(self, chrom, hapl_type):
        bks = filter(lambda item: item.chrom == chrom and item.hapl_type ==
                     hapl_type, self.breaks_list)
        if len(bks) > 0:
            return True
        else:
            return False

    def _generateHapStr(self, chrom, hap_i, ref):
        hapStr = ""

        bks = filter(lambda item: item.chrom == chrom and item.hapl_type ==
                     hap_i, self.breaks_list)
        bks = sorted(bks,
                     key=lambda item: (item.position, item.paired_position))
        posis = [0] + [bks[i].position for i in range(len(bks))] \
            + [len(ref[chrom][hap_i])]
        refSegs = [ref[chrom][hap_i][posis[i]:posis[i+1]] for i in
                   range(len(posis) - 1)]

        delStart = False
        for i in range(len(refSegs)):
            if not delStart:
                hapStr = hapStr + refSegs[i]
            if i < len(bks):
                if bks[i].name == "duplication":
                    hapStr = hapStr + bks[i].insertStr
                elif bks[i].name == "deletion":
                    if bks[i].position < bks[i].paired_position:
                        if delStart:
                            print "Error"
                        delStart = True
                    else:
                        if not delStart:
                            print "Error"
                        delStart = False

        return hapStr

    def _getRandomPois(self, chrom, avail_position):
        pois = -1
        while True:
            # 从候选位置中抽取插入位置，位置pois - 1
            pois = avail_position.sample1pois(chrom)
            print "pois sampled :{}".format(pois)
            if avail_position.isOverlaped(chrom, pois-1, pois):
                avail_position.takePois(chrom, pois-1, pois)
                break
            else:
                continue

        return pois

    def _generateDupBPs(self, variant_dict, avail_position, ref, ploidy_status):
        """
        生成duplication copy
        """
        # dic value: 位置、变异长度、拷贝数、基因型、名称
        # 注意生成多重ploidy的情况

        for chrom in variant_dict.keys():
            psc = Counter(ploidy_status[chrom])
            for vartItem in variant_dict[chrom]:
                if vartItem[-2] == "NONE":
                    continue
                else:
                    genoCounter = Counter(vartItem[-2])

                    hap_set_genoc = set(genoCounter.keys())
                    hap_set_plosc = set(psc.keys())
                    assert hap_set_genoc <= hap_set_plosc

                    for genoHap in genoCounter.keys():
                        for copy_i in range(genoCounter[genoHap] -
                                            psc[genoHap]):
                            bp = Break()
                            bp.chrom = chrom
                            bp.hapl_type = random.sample(['P', 'M'], 1)[0]
                            bp.hapl_index = random.sample(
                                range(psc[genoHap]), 1)[0]

                            bp.position = self._getRandomPois(
                                chrom, avail_position)

                            bp.name = "duplication"
                            if genoHap == 'P':
                                bp.insertStr = ref[chrom]['P'][hap][
                                    vartItem[0]: vartItem[0] + vartItem[1]]
                            elif genoHap == 'M':
                                bp.insertStr = ref[chrom]['M'][hap][
                                    vartItem[0]: vartItem[0] + vartItem[1]]
                            else:
                                print "Error"
                            self.breaks_list.append(bp)

    def _generateDelBPs(self, variant_dict, ploidy_status):
        """ 生成deletion的breakpoints
        :returns: TODO

        """
        # dic value: 位置、变异长度、拷贝数、基因型、名称

        for chrom in variant_dict.keys():
            # 获取当前染色体的ploidy状态
            psc = Counter(ploidy_status[chrom])

            for vartItem in variant_dict[chrom]:
                genotype = vartItem[-2]
                if  genotype == "NONE":
                    # 每一个单体都要形成一对del bp
                    for hapl_type in psc.keys():
                        for hapl_index in range(psc[hapl_type]):
                            bp_start = self._deletionBP(
                                chrom, hapl_type, hapl_index,
                                vartItem[0], vartItem[0]+vartItem[1])
                            bp_end = self._deletionBP(
                                chrom, hapl_type, hapl_index,
                                vartItem[0]+vartItem[1], vartItem[0])
                            self.breaks_list = self.breaks_list + \
                                [bp_start, bp_end]
                else:
                    genoCounter = Counter(genotype)
                    # 只要出现与ploidy状态，减少的位置就要写入deletion

                    # 出现纯合时候, No key
                    if len(genoCounter.keys()) < len(psc.keys()):
                        for hapl_type in set(psc.keys)-set(genoCounter.keys()):
                            for hapl_index in range(psc[hapl_type]):
                                bp_start = self._deletionBP(
                                    chrom, hapl_type, hapl_index,
                                    vartItem[0], vartItem[0]+vartItem[1])
                                bp_end = self._deletionBP(
                                    chrom, hapl_type, hapl_index,
                                    vartItem[0]+vartItem[1], vartItem[0])
                                self.breaks_list = self.breaks_list + \
                                    [bp_start, bp_end]

                    for hapl_type in genoCounter.keys():
                        del_cp_num = psc[hapl_type] - genoCounter[hapl_type]
                        if del_cp_num <= 0:
                            continue
                        else:
                            del_cp_indexes = ramdom.sample(
                                range(psc[hapl_type]), del_cp_num)
                            for hapl_index in del_cp_indexes:
                                bp_start = self._deletionBP(
                                    chrom, hapl_type, hapl_index,
                                    vartItem[0], vartItem[0]+vartItem[1])
                                bp_end = self._deletionBP(
                                    chrom, hapl_type, hapl_index,
                                    vartItem[0]+vartItem[1], vartItem[0])
                                self.breaks_list = self.breaks_list + \
                                    [bp_start, bp_end]


    def _deletionBP(
            self,
            chrom,
            hapl_type,
            hapl_index,
            position,
            paired_position):
        """TODO: Docstring for _generateBPdeletion.

        :arg1: TODO
        :returns: TODO

        """
        bp = Break()
        bp.name = "deletion"

        bp.chrom = chrom
        bp.hapl_type = hapl_type
        bp.hapl_index = hapl_index
        bp.position = position
        bp.paired_position = paired_position

        return bp
