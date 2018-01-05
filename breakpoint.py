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
from utils.utils import rand_DNA


class Break:

    def __init__(self):
        self.chrom = ""
        self.hapl_type = -1
        # it seems , no need
        self.hapl_idx = -1
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

        # 此处为了进行ploidy操作，结构为 {chr1:{'P':[[bp1],[bp2]], 'M':[[bp3]]}}
        self.breaks_list = {}

    def add_ploidy(self, chrom, hapl_type, hapl_idxes):
        hic = Counter(hapl_idxes)
        ploidies = []
        for hi in hic.keys():
            hapl_idx = int(hi)
            number = hic[hi]
            if self._has(chrom, hapl_type, hapl_idx):
                ploidies = ploidies +\
                    number * [self.breaks_list[chrom][hapl_type][hapl_idx]]

        if self._has(chrom, hapl_type, hapl_idx):
            self.breaks_list[chrom][hapl_type] = self.breaks_list[chrom][
                hapl_type] + ploidies

    def _has(self, chrom, hapl_type, hapl_idx):
        if chrom in self.breaks_list.keys():
            if hapl_type in self.breaks_list[chrom].keys():
                if hapl_idx < len(self.breaks_list[chrom][hapl_type]):
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

        # self.breakpoints.delete_ploidy(chrom, hapl, hapl_idxes)
    def delete_ploidy(self, chrom, hapl, hapl_idxes):
        self.breaks_list[chrom][hapl] = [
            self.breaks_list[chrom][hapl][i]
            for i in range(len(self.breaks_list[chrom][hapl]))
            if i not in hapl_idxes]

    def generateBPs(self, variant_positions, avail_position, noDEL_position,
                    ref, ploidy_status):

        for chrom in variant_positions.svp_dict.keys():
            # CNV
            vps_CNV = filter(lambda item: item.sv_type == "CNV",
                             variant_positions.svp_dict[chrom])
            self._generateCDupBPs(chrom, vps_CNV, avail_position,
                                  noDEL_position, ref,
                                  ploidy_status)
            self._generateCDelBPs(chrom, vps_CNV, noDEL_position, ploidy_status)

            # INSERTION
            vp_INSERTION = filter(
                lambda item: item.sv_type == "INSERTION",
                variant_positions.svp_dict[chrom])
            self._generateInsBPs(chrom, vp_INSERTION,
                                 avail_position, ploidy_status)

            # DELETION
            vp_DELETION = filter(
                lambda item: item.sv_type == "DELETION",
                variant_positions.svp_dict[chrom])
            self._generateDelsBPs(chrom, vp_DELETION, ploidy_status)

            # COMPLEXINDEL
            vp_COMPLEXINDEL = filter(
                lambda item: item.sv_type == "COMPLEXINDEL",
                variant_positions.svp_dict[chrom])
            self._generateComplexindelsBPs(chrom, vp_COMPLEXINDEL, ploidy_status)

            vp_INVERTION = filter(
                lambda item: item.sv_type == "INVERTION",
                variant_positions.svp_dict[chrom])
            self._generateInvsBPs(chrom, vp_INVERTION,
                                  ploidy_status,
                                  ref)

            vp_TANDEMDUP = filter(
                lambda item: item.sv_type == "TANDEMDUP",
                variant_positions.svp_dict[chrom])
            self._generateTandemdupBPs(chrom, vp_TANDEMDUP, ploidy_status, ref)

            vp_TRANSLOCATION = filter(
                lambda item: item.sv_type == "TRANSLOCATION",
                variant_positions.svp_dict[chrom])
            self._generateTransBPs(chrom, vp_TRANSLOCATION, avail_position,
                                   ploidy_status, ref)

    def _generateTransBPs(self, chrom, vps, avail_position, ploidy_status, ref):
        for vp in vps:
            position = vp.position
            length = vp.sv.length
            htf = vp.sv.hapl_type_from
            htt = vp.sv.hapl_type_to
            hif = vp.sv.hapl_idx_from
            hit = vp.sv.hapl_idx_to

            bp_start = self._pairedBP("DELETION", chrom, htf, hif,
                                      position, position+length)
            bp_end = self._pairedBP("DELETION", chrom, htf, hif,
                                    position+length, position)

            self._breakAppend(ploidy_status, chrom, htf, hif, bp_start)
            self._breakAppend(ploidy_status, chrom, htf, hif, bp_end)

            bp = Break()
            bp.chrom = vp.sv.chrom_to
            bp.hapl_type = htt
            bp.hapl_idx = hit
            bp.position = self._getRandomPosi(bp.chrom, avail_position)
            bp.name = "INSERTION"
            bp.insertStr = ref[chrom][htf][hif][position-1: position+length]
            self._breakAppend(ploidy_status, chrom, htt, hit, bp)

    def _generateInvsBPs(self, chrom, vps, ploidy_status, ref):
        for vp in vps:
            hapl_type = vp.sv.hapl_type
            hapl_idx = vp.sv.hapl_idx

            insertStr =\
                ref[chrom][hapl_type][hapl_idx][
                    vp.position: vp.position+vp.sv.length][::-1]

            bp_start = self._pairedBP(
                "INVERTION",
                chrom, hapl_type, hapl_idx,
                vp.position,
                vp.position+vp.sv.length,
                insertStr)
            bp_end = self._pairedBP("INVERTION",
                                    chrom, hapl_type, hapl_idx,
                                    vp.position+vp.sv.length,
                                    vp.position,
                                    insertStr)

            self._breakAppend(ploidy_status, chrom,
                              hapl_type, hapl_idx, bp_start)
            self._breakAppend(ploidy_status, chrom,
                              hapl_type, hapl_idx, bp_end)

    def _generateTandemdupBPs(self, chrom, vps, ploidy_status, ref):
        for vp in vps:
            hapl_type = vp.sv.hapl_type
            hapl_idx = vp.sv.hapl_idx
            times = vp.sv.times

            insertStr =\
                times * ref[chrom][hapl_type][hapl_idx][
                    vp.position: vp.position+vp.sv.length]

            bp = Break()
            bp.chrom = chrom
            bp.hapl_type = hapl_type
            bp.hapl_idx = hapl_idx
            bp.position = vp.position+vp.sv.length
            bp.name = "INSERTION"
            bp.insertStr = insertStr

            self._breakAppend(ploidy_status, chrom, hapl_type, hapl_idx, bp)

    def _generateDelsBPs(self, chrom, vps, ploidy_status):
        for vp in vps:

            hapl_type = vp.sv.hapl_type
            hapl_idx = vp.sv.hapl_idx

            bp_start = self._pairedBP(
                "DELETION",
                chrom, hapl_type, hapl_idx,
                vp.position,
                vp.position+vp.sv.length)
            bp_end = self._pairedBP("DELETION",
                                    chrom, hapl_type, hapl_idx,
                                    vp.position+vp.sv.length,
                                    vp.position)

            self._breakAppend(ploidy_status, chrom,
                              hapl_type, hapl_idx, bp_start)
            self._breakAppend(ploidy_status, chrom,
                              hapl_type, hapl_idx, bp_end)

    def _generateComplexindelsBPs(self, chrom, vps, ploidy_status):
        for vp in vps:
            # self.length = -1
            # self.hapl_type = ""
            # self.hapl_idx = -1
            hapl_type = vp.sv.hapl_type
            hapl_idx = vp.sv.hapl_idx
            insertStr = rand_DNA(vp.sv.length2)

            bp_start = self._pairedBP(
                "COMPLEXINDEL",
                chrom, hapl_type, hapl_idx,
                vp.position,
                vp.position+vp.sv.length1, insertStr)
            bp_end = self._pairedBP("COMPLEXINDEL",
                                    chrom, hapl_type, hapl_idx,
                                    vp.position+vp.sv.length1,
                                    vp.position, insertStr)

            self._breakAppend(ploidy_status, chrom,
                              hapl_type, hapl_idx, bp_start)
            self._breakAppend(ploidy_status, chrom,
                              hapl_type, hapl_idx, bp_end)

    def _generateInsBPs(self, chrom, vps, avail_position, ploidy_status):
        for vp in vps:
            # self.length = -1
            # self.hapl_type = ""
            # self.hapl_idx = -1
            hapl_type = vp.sv.hapl_type
            hapl_idx = vp.sv.hapl_idx
            length = vp.sv.length

            position = vp.position

            insertStr = rand_DNA(length)

            bp = Break()
            bp.chrom = chrom
            bp.hapl_type = hapl_type
            bp.hapl_idx = hapl_idx
            bp.position = position
            bp.name = "INSERTION"
            bp.insertStr = insertStr
            self._breakAppend(ploidy_status, chrom,
                              hapl_type, hapl_idx, bp)

    def generateFa(self, ref):
        newFa = copy.deepcopy(ref)
        for chrom in ref.keys():
            if not self._bpsHasChrom(chrom):
                continue
            else:
                for hapl_type in ref[chrom].keys():
                    if not self._bpsHasChromHap(chrom, hapl_type):
                        continue
                    else:
                        for hapl_idx in\
                                range(len(self.breaks_list[chrom][hapl_type])):
                            newFa[chrom][hapl_type][hapl_idx] =\
                                self._generateHapStr(chrom, hapl_type, hapl_idx,
                                                     ref)
        return newFa

    def _bpsHasChrom(self, chrom):
        return chrom in self.breaks_list.keys()

    def _bpsHasChromHap(self, chrom, hapl_type):
        return self._bpsHasChrom(chrom) and hapl_type in\
            self.breaks_list[chrom].keys()

    def _generateHapStr(self, chrom, hapl_type, hapl_idx, ref):
        hapStr = ""

        bks = self.breaks_list[chrom][hapl_type][hapl_idx]
        bks = sorted(bks,
                     key=lambda item: (item.position, item.paired_position))
        posis = [0] + [bks[i].position for i in range(len(bks))] \
            + [len(ref[chrom][hapl_type][hapl_idx])]
        refSegs = [ref[chrom][hapl_type][hapl_idx][posis[i]:posis[i+1]] for i in
                   range(len(posis) - 1)]

        delStart = False
        invStart = False
        for i in range(len(refSegs)):
            if not delStart:
                hapStr = hapStr + refSegs[i]
            if i < len(bks):
                if bks[i].name == "DUPLICATION" or\
                        bks[i].name == "INSERTION" or\
                        bks[i].name == "TRANSLOCATION":
                    hapStr = hapStr + bks[i].insertStr
                if bks[i].name == "DELETION":
                    if bks[i].position < bks[i].paired_position:
                        if delStart:
                            print "Error"
                        delStart = True
                    else:
                        if not delStart:
                            print "Error"
                        delStart = False
                if bks[i].name == "INVERTION":
                    if bks[i].position < bks[i].paired_position:
                        if invStart:
                            print "Error"
                        invStart = True
                    else:
                        hapStr = hapStr + bks[i].insertStr
                        if not invStart:
                            print "Error"
                        invStart = False

        return hapStr

    def _getRandomPosi(self, chrom, avail_position):
        pois = -1
        while True:
            # 从候选位置中抽取插入位置，位置pois - 1
            pois = avail_position.sample1posi(chrom)
            print "pois sampled :{}".format(pois)
            if avail_position.isOverlaped(chrom, pois, pois+1):
                avail_position.takePosi(chrom, pois, pois+1)
                break
            else:
                continue

        return pois

    def _getRandomInvalPosi(self, chrom, length, avail_position):
        pois = -1
        while True:
            # 从候选位置中抽取插入位置，位置pois
            pois = avail_position.sample1posi(chrom)
            print "pois sampled :{}".format(pois)
            if avail_position.isOverlaped(chrom, pois, pois+length):
                avail_position.takePosi(chrom, pois, pois+length)
                break
            else:
                continue

        return pois

    def _generateCDupBPs(self, chrom, vps, avail_position, noDEL_position,
                         ref, ploidy_status):
        """
        生成duplication copy
        """
        # 此时的ref是经过了ploidy和snv之后的，以便于生成cnv

        # dic value: 位置、变异长度、拷贝数、基因型、名称
        # 注意生成多重ploidy的情况

        psc = Counter(ploidy_status[chrom])
        for vp in vps:
            genotype = vp.sv.genotype
            if genotype == "NONE":
                continue
            else:
                noDEL_position.addRange(chrom, vp.position,
                                        vp.position+vp.sv.length)

                genoCounter = Counter(genotype)

                hap_set_genoc = set(genoCounter.keys())
                hap_set_plosc = set(psc.keys())

                for genoHap in hap_set_genoc & hap_set_plosc:
                    for copy_i in range(genoCounter[genoHap] -
                                        psc[genoHap]):
                        bp = Break()
                        bp.chrom = chrom
                        bp.hapl_type = random.sample(['P', 'M'], 1)[0]
                        bp.hapl_idx = random.sample(
                            range(len(ref[chrom][bp.hapl_type])), 1)[0]

                        bp.position = self._getRandomPosi(
                            chrom, avail_position)

                        bp.name = "DUPLICATION"

                        hapi = random.sample(range(psc[genoHap]), 1)[0]

                        bp.insertStr = ref[chrom][genoHap][hapi][
                            vp.position:vp.position+vp.sv.length]

                        self._breakAppend(ploidy_status, chrom,
                                          bp.hapl_type, bp.hapl_idx, bp)

    def _breakAppend(self, ploidy_status, chrom, hapl_type, hapl_idx, bp):
        # 此处为了进行ploidy操作，结构为 {chr1:{'P':[[bp1],[bp2]], 'M':[[bp3]]}}
        psc = Counter(ploidy_status[chrom])
        if chrom not in self.breaks_list.keys():
            self.breaks_list[chrom] = {}
            for hap_key in psc.keys():
                self.breaks_list[chrom][hap_key] = []
                for i in range(psc[hap_key]):
                    self.breaks_list[chrom][hap_key].append([])
        else:
            if hapl_type not in self.breaks_list[chrom].keys():
                self.breaks_list[chrom][hapl_type] = []
                for i in range(psc[hapl_type]):
                    self.breaks_list[chrom][hapl_type].append([])
            else:
                if hapl_idx not in range(
                        len(self.breaks_list[chrom][hapl_type])):
                    for i in range(hapl_idx -
                                   len(self.breaks_list[chrom][hapl_type]) + 1):
                        self.breaks_list[chrom][hapl_type].append([])

        self.breaks_list[chrom][hapl_type][hapl_idx].append(bp)

    def _generateCDelBPs(self, chrom, vps, noDEL_position, ploidy_status):
        """ 生成deletion的breakpoints
        :returns: TODO

        """
        # dic value: 位置、变异长度、拷贝数、基因型、名称

        psc = Counter(ploidy_status[chrom])

        for vp in vps:
            genotype = vp.sv.genotype
            if genotype == "NONE":
                # 每一个单体都要形成一对del bp
                for hapl_type in psc.keys():
                    for hapl_idx in range(psc[hapl_type]):
                        bp_start = self._pairedBP("DELETION",
                                                  chrom, hapl_type, hapl_idx,
                                                  vp.position,
                                                  vp.position+vp.sv.length)
                        bp_end = self._pairedBP("DELETION",
                                                chrom, hapl_type, hapl_idx,
                                                vp.position+vp.sv.length,
                                                vp.position)

                        self._breakAppend(ploidy_status, chrom,
                                          hapl_type, hapl_idx, bp_start)
                        self._breakAppend(ploidy_status, chrom,
                                          hapl_type, hapl_idx, bp_end)
            else:
                genoCounter = Counter(genotype)
                # 只要出现与ploidy状态，减少的位置就要写入deletion

                # 出现纯合时候, No key
                if len(genoCounter.keys()) < len(psc.keys()):
                    for hapl_type in set(psc.keys())-set(genoCounter.keys()):
                        for hapl_idx in range(psc[hapl_type]):
                            bp_start = self._pairedBP("DELETION",
                                                      chrom, hapl_type,
                                                      hapl_idx,
                                                      vp.position,
                                                      vp.position+vp.sv.length)
                            bp_end = self._pairedBP("DELETION",
                                                    chrom, hapl_type, hapl_idx,
                                                    vp.position+vp.sv.length,
                                                    vp.position)
                            self._breakAppend(ploidy_status, chrom,
                                              hapl_type, hapl_idx, bp_start)
                            self._breakAppend(ploidy_status, chrom,
                                              hapl_type, hapl_idx, bp_end)

                for hapl_type in genoCounter.keys():
                    del_cp_num = psc[hapl_type] - genoCounter[hapl_type]
                    if del_cp_num <= 0:
                        continue
                    else:
                        del_cp_indexes = random.sample(
                            range(psc[hapl_type]), del_cp_num)
                        for hapl_idx in del_cp_indexes:
                            bp_start = self._pairedBP("DELETION",
                                                      chrom, hapl_type,
                                                      hapl_idx,
                                                      vp.position,
                                                      vp.position+vp.sv.length)
                            bp_end = self._pairedBP("DELETION",
                                                    chrom, hapl_type, hapl_idx,
                                                    vp.position+vp.sv.length,
                                                    vp.position)
                            self._breakAppend(ploidy_status, chrom,
                                              hapl_type, hapl_idx, bp_start)
                            self._breakAppend(ploidy_status, chrom,
                                              hapl_type, hapl_idx, bp_end)

    def _pairedBP(
            self,
            name,
            chrom,
            hapl_type,
            hapl_idx,
            position,
            paired_position, insertStr=None):
        """TODO: Docstring for _generateBPdeletion.

        :arg1: TODO
        :returns: TODO

        """
        bp = Break()
        bp.name = name
        bp.chrom = chrom
        bp.hapl_type = hapl_type
        bp.hapl_idx = hapl_idx
        bp.position = position
        bp.paired_position = paired_position
        bp.insertStr = insertStr

        return bp
