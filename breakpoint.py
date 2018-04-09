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
        self.haplType = -1
        # it seems , no need
        self.haplIdx = -1
        self.position = -1
        self.name = ""
        self.insertStr = None

        # deletion的paired position
        # l
        self.pairedPosition = -1


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
        self.breaksList = {}

    def add_ploidy(self, chrom, haplType, hapl_idxes):
        hic = Counter(hapl_idxes)
        ploidies = []
        for hi in hic.keys():
            haplIdx = int(hi)
            number = hic[hi]
            if self._has(chrom, haplType, haplIdx):
                ploidies = ploidies +\
                    number * [self.breaksList[chrom][haplType][haplIdx]]

        if self._has(chrom, haplType, haplIdx):
            self.breaksList[chrom][haplType] = self.breaksList[chrom][
                haplType] + ploidies

    def _has(self, chrom, haplType, haplIdx):
        if chrom in self.breaksList.keys():
            if haplType in self.breaksList[chrom].keys():
                if haplIdx < len(self.breaksList[chrom][haplType]):
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

        # self.breakpoints.delete_ploidy(chrom, hapl, hapl_idxes)
    def delete_ploidy(self, chrom, hapl, hapl_idxes):
        self.breaksList[chrom][hapl] = [
            self.breaksList[chrom][hapl][i]
            for i in range(len(self.breaksList[chrom][hapl]))
            if i not in hapl_idxes]

    def generateBPs(self, variantPositions, availPosition, noDELPosition,
                    ref, ploidyStatus):

        for chrom in variantPositions.svpDict.keys():
            # CNV
            vpsCNV = filter(lambda item: item.sv_type == "CNV",
                             variantPositions.svpDict[chrom])
            self._generateCDupBPs(chrom, vpsCNV, availPosition,
                                  noDELPosition, ref,
                                  ploidyStatus)
            self._generateCDelBPs(chrom, vpsCNV, noDELPosition, ploidyStatus)

            # INSERTION
            vpINSERTION = filter(
                lambda item: item.sv_type == "INSERTION",
                variantPositions.svpDict[chrom])
            self._generateInsBPs(chrom, vpINSERTION,
                                 availPosition, ploidyStatus)

            # DELETION
            vpDELETION = filter(
                lambda item: item.sv_type == "DELETION",
                variantPositions.svpDict[chrom])
            self._generateDelsBPs(chrom, vpDELETION, ploidyStatus)

            # COMPLEXINDEL
            vpCOMPLEXINDEL = filter(
                lambda item: item.sv_type == "COMPLEXINDEL",
                variantPositions.svpDict[chrom])
            self._generateComplexindelsBPs(chrom, vpCOMPLEXINDEL, ploidyStatus)

            vpINVERTION = filter(
                lambda item: item.sv_type == "INVERTION",
                variantPositions.svpDict[chrom])
            self._generateInvsBPs(chrom, vpINVERTION,
                                  ploidyStatus,
                                  ref)

            vpTANDEMDUP = filter(
                lambda item: item.sv_type == "TANDEMDUP",
                variantPositions.svpDict[chrom])
            self._generateTandemdupBPs(chrom, vpTANDEMDUP, ploidyStatus, ref)

            vpTRANSLOCATION = filter(
                lambda item: item.sv_type == "TRANSLOCATION",
                variantPositions.svpDict[chrom])
            self._generateTransBPs(chrom, vpTRANSLOCATION, availPosition,
                                   ploidyStatus, ref)

    def _generateTransBPs(self, chrom, vps, availPosition, ploidyStatus, ref):
        for vp in vps:
            position = vp.position
            length = vp.sv.length
            htf = vp.sv.haplTypeFrom
            htt = vp.sv.haplTypeTo
            hif = vp.sv.haplIdxFrom
            hit = vp.sv.haplIdxTo

            bpStart = self._pairedBP("DELETION", chrom, htf, hif,
                                      position, position+length)
            bpEnd = self._pairedBP("DELETION", chrom, htf, hif,
                                    position+length, position)

            self._breakAppend(ploidyStatus, chrom, htf, hif, bpStart)
            self._breakAppend(ploidyStatus, chrom, htf, hif, bpEnd)

            bp = Break()
            bp.chrom = vp.sv.chrom_to
            bp.haplType = htt
            bp.haplIdx = hit
            bp.position = self._getRandomPosi(bp.chrom, availPosition)
            bp.name = "INSERTION"
            bp.insertStr = ref[chrom][htf][hif][position-1: position+length]
            self._breakAppend(ploidyStatus, chrom, htt, hit, bp)

    def _generateInvsBPs(self, chrom, vps, ploidyStatus, ref):
        for vp in vps:
            haplType = vp.sv.haplType
            haplIdx = vp.sv.haplIdx

            insertStr =\
                ref[chrom][haplType][haplIdx][
                    vp.position: vp.position+vp.sv.length][::-1]

            bpStart = self._pairedBP(
                "INVERTION",
                chrom, haplType, haplIdx,
                vp.position,
                vp.position+vp.sv.length,
                insertStr)
            bpEnd = self._pairedBP("INVERTION",
                                    chrom, haplType, haplIdx,
                                    vp.position+vp.sv.length,
                                    vp.position,
                                    insertStr)

            self._breakAppend(ploidyStatus, chrom,
                              haplType, haplIdx, bpStart)
            self._breakAppend(ploidyStatus, chrom,
                              haplType, haplIdx, bpEnd)

    def _generateTandemdupBPs(self, chrom, vps, ploidyStatus, ref):
        for vp in vps:
            haplType = vp.sv.haplType
            haplIdx = vp.sv.haplIdx
            times = vp.sv.times

            insertStr =\
                times * ref[chrom][haplType][haplIdx][
                    vp.position: vp.position+vp.sv.length]

            bp = Break()
            bp.chrom = chrom
            bp.haplType = haplType
            bp.haplIdx = haplIdx
            bp.position = vp.position+vp.sv.length
            bp.name = "INSERTION"
            bp.insertStr = insertStr

            self._breakAppend(ploidyStatus, chrom, haplType, haplIdx, bp)

    def _generateDelsBPs(self, chrom, vps, ploidyStatus):
        for vp in vps:

            haplType = vp.sv.haplType
            haplIdx = vp.sv.haplIdx

            bpStart = self._pairedBP(
                "DELETION",
                chrom, haplType, haplIdx,
                vp.position,
                vp.position+vp.sv.length)
            bpEnd = self._pairedBP("DELETION",
                                    chrom, haplType, haplIdx,
                                    vp.position+vp.sv.length,
                                    vp.position)

            self._breakAppend(ploidyStatus, chrom,
                              haplType, haplIdx, bpStart)
            self._breakAppend(ploidyStatus, chrom,
                              haplType, haplIdx, bpEnd)

    def _generateComplexindelsBPs(self, chrom, vps, ploidyStatus):
        for vp in vps:
            # self.length = -1
            # self.haplType = ""
            # self.haplIdx = -1
            haplType = vp.sv.haplType
            haplIdx = vp.sv.haplIdx
            insertStr = rand_DNA(vp.sv.length2)

            bpStart = self._pairedBP(
                "COMPLEXINDEL",
                chrom, haplType, haplIdx,
                vp.position,
                vp.position+vp.sv.length1, insertStr)
            bpEnd = self._pairedBP("COMPLEXINDEL",
                                    chrom, haplType, haplIdx,
                                    vp.position+vp.sv.length1,
                                    vp.position, insertStr)

            self._breakAppend(ploidyStatus, chrom,
                              haplType, haplIdx, bpStart)
            self._breakAppend(ploidyStatus, chrom,
                              haplType, haplIdx, bpEnd)

    def _generateInsBPs(self, chrom, vps, availPosition, ploidyStatus):
        for vp in vps:
            # self.length = -1
            # self.haplType = ""
            # self.haplIdx = -1
            haplType = vp.sv.haplType
            haplIdx = vp.sv.haplIdx
            length = vp.sv.length

            position = vp.position

            insertStr = rand_DNA(length)

            bp = Break()
            bp.chrom = chrom
            bp.haplType = haplType
            bp.haplIdx = haplIdx
            bp.position = position
            bp.name = "INSERTION"
            bp.insertStr = insertStr
            self._breakAppend(ploidyStatus, chrom,
                              haplType, haplIdx, bp)

    def generateFa(self, ref):
        newFa = copy.deepcopy(ref)
        for chrom in ref.keys():
            if not self._bpsHasChrom(chrom):
                continue
            else:
                for haplType in ref[chrom].keys():
                    if not self._bpsHasChromHap(chrom, haplType):
                        continue
                    else:
                        for haplIdx in\
                                range(len(self.breaksList[chrom][haplType])):
                            newFa[chrom][haplType][haplIdx] =\
                                self._generateHapStr(chrom, haplType, haplIdx,
                                                     ref)
        return newFa

    def _bpsHasChrom(self, chrom):
        return chrom in self.breaksList.keys()

    def _bpsHasChromHap(self, chrom, haplType):
        return self._bpsHasChrom(chrom) and haplType in\
            self.breaksList[chrom].keys()

    def _generateHapStr(self, chrom, haplType, haplIdx, ref):
        hapStr = ""

        bks = self.breaksList[chrom][haplType][haplIdx]
        bks = sorted(bks,
                     key=lambda item: (item.position, item.pairedPosition))
        posis = [0] + [bks[i].position for i in range(len(bks))] \
            + [len(ref[chrom][haplType][haplIdx])]
        refSegs = [ref[chrom][haplType][haplIdx][posis[i]:posis[i+1]] for i in
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
                    if bks[i].position < bks[i].pairedPosition:
                        if delStart:
                            print "Error"
                        delStart = True
                    else:
                        if not delStart:
                            print "Error"
                        delStart = False
                if bks[i].name == "INVERTION":
                    if bks[i].position < bks[i].pairedPosition:
                        if invStart:
                            print "Error"
                        invStart = True
                    else:
                        hapStr = hapStr + bks[i].insertStr
                        if not invStart:
                            print "Error"
                        invStart = False

        return hapStr

    def _getRandomPosi(self, chrom, availPosition):
        pois = -1
        while True:
            # 从候选位置中抽取插入位置，位置pois - 1
            pois = availPosition.sample1posi(chrom)
            print "pois sampled :{}".format(pois)
            if availPosition.isOverlaped(chrom, pois, pois+1):
                availPosition.takePosi(chrom, pois, pois+1)
                break
            else:
                continue

        return pois

    def _getRandomInvalPosi(self, chrom, length, availPosition):
        pois = -1
        while True:
            # 从候选位置中抽取插入位置，位置pois
            pois = availPosition.sample1posi(chrom)
            print "pois sampled :{}".format(pois)
            if availPosition.isOverlaped(chrom, pois, pois+length):
                availPosition.takePosi(chrom, pois, pois+length)
                break
            else:
                continue

        return pois

    def _generateCDupBPs(self, chrom, vps, availPosition, noDELPosition,
                         ref, ploidyStatus):
        """
        生成duplication copy
        """
        # 此时的ref是经过了ploidy和snv之后的，以便于生成cnv

        # dic value: 位置、变异长度、拷贝数、基因型、名称
        # 注意生成多重ploidy的情况

        psc = Counter(ploidyStatus[chrom])
        for vp in vps:
            genotype = vp.sv.genotype
            if genotype == "NONE":
                continue
            else:
                noDELPosition.addRange(chrom, vp.position,
                                        vp.position+vp.sv.length)

                genoCounter = Counter(genotype)

                hapSetGenoc = set(genoCounter.keys())
                hapSetPlosc = set(psc.keys())

                for genoHap in hapSetGenoc & hapSetPlosc:
                    for copy_i in range(genoCounter[genoHap] -
                                        psc[genoHap]):
                        bp = Break()
                        bp.chrom = chrom
                        bp.haplType = random.sample(['P', 'M'], 1)[0]
                        bp.haplIdx = random.sample(
                            range(len(ref[chrom][bp.haplType])), 1)[0]

                        bp.position = self._getRandomPosi(
                            chrom, availPosition)

                        bp.name = "DUPLICATION"

                        hapi = random.sample(range(psc[genoHap]), 1)[0]

                        bp.insertStr = ref[chrom][genoHap][hapi][
                            vp.position:vp.position+vp.sv.length]

                        self._breakAppend(ploidyStatus, chrom,
                                          bp.haplType, bp.haplIdx, bp)

    def _breakAppend(self, ploidyStatus, chrom, haplType, haplIdx, bp):
        # 此处为了进行ploidy操作，结构为 {chr1:{'P':[[bp1],[bp2]], 'M':[[bp3]]}}
        psc = Counter(ploidyStatus[chrom])
        if chrom not in self.breaksList.keys():
            self.breaksList[chrom] = {}
            for hap_key in psc.keys():
                self.breaksList[chrom][hap_key] = []
                for i in range(psc[hap_key]):
                    self.breaksList[chrom][hap_key].append([])
        else:
            if haplType not in self.breaksList[chrom].keys():
                self.breaksList[chrom][haplType] = []
                for i in range(psc[haplType]):
                    self.breaksList[chrom][haplType].append([])
            else:
                if haplIdx not in range(
                        len(self.breaksList[chrom][haplType])):
                    for i in range(haplIdx -
                                   len(self.breaksList[chrom][haplType]) + 1):
                        self.breaksList[chrom][haplType].append([])

        self.breaksList[chrom][haplType][haplIdx].append(bp)

    def _generateCDelBPs(self, chrom, vps, noDELPosition, ploidyStatus):
        """ 生成deletion的breakpoints
        :returns: TODO

        """
        # dic value: 位置、变异长度、拷贝数、基因型、名称

        psc = Counter(ploidyStatus[chrom])

        for vp in vps:
            genotype = vp.sv.genotype
            if genotype == "NONE":
                # 每一个单体都要形成一对del bp
                for haplType in psc.keys():
                    for haplIdx in range(psc[haplType]):
                        bpStart = self._pairedBP("DELETION",
                                                  chrom, haplType, haplIdx,
                                                  vp.position,
                                                  vp.position+vp.sv.length)
                        bpEnd = self._pairedBP("DELETION",
                                                chrom, haplType, haplIdx,
                                                vp.position+vp.sv.length,
                                                vp.position)

                        self._breakAppend(ploidyStatus, chrom,
                                          haplType, haplIdx, bpStart)
                        self._breakAppend(ploidyStatus, chrom,
                                          haplType, haplIdx, bpEnd)
            else:
                genoCounter = Counter(genotype)
                # 只要出现与ploidy状态，减少的位置就要写入deletion

                # 出现纯合时候, No key
                if len(genoCounter.keys()) < len(psc.keys()):
                    for haplType in set(psc.keys())-set(genoCounter.keys()):
                        for haplIdx in range(psc[haplType]):
                            bpStart = self._pairedBP("DELETION",
                                                      chrom, haplType,
                                                      haplIdx,
                                                      vp.position,
                                                      vp.position+vp.sv.length)
                            bpEnd = self._pairedBP("DELETION",
                                                    chrom, haplType, haplIdx,
                                                    vp.position+vp.sv.length,
                                                    vp.position)
                            self._breakAppend(ploidyStatus, chrom,
                                              haplType, haplIdx, bpStart)
                            self._breakAppend(ploidyStatus, chrom,
                                              haplType, haplIdx, bpEnd)

                for haplType in genoCounter.keys():
                    delCpNum = psc[haplType] - genoCounter[haplType]
                    if delCpNum <= 0:
                        continue
                    else:
                        delCpIndexes = random.sample(
                            range(psc[haplType]), delCpNum)
                        for haplIdx in delCpIndexes:
                            bpStart = self._pairedBP("DELETION",
                                                      chrom, haplType,
                                                      haplIdx,
                                                      vp.position,
                                                      vp.position+vp.sv.length)
                            bpEnd = self._pairedBP("DELETION",
                                                    chrom, haplType, haplIdx,
                                                    vp.position+vp.sv.length,
                                                    vp.position)
                            self._breakAppend(ploidyStatus, chrom,
                                              haplType, haplIdx, bpStart)
                            self._breakAppend(ploidyStatus, chrom,
                                              haplType, haplIdx, bpEnd)

    def _pairedBP(
            self,
            name,
            chrom,
            haplType,
            haplIdx,
            position,
            pairedPosition, insertStr=None):
        """TODO: Docstring for _generateBPdeletion.

        :arg1: TODO
        :returns: TODO

        """
        bp = Break()
        bp.name = name
        bp.chrom = chrom
        bp.haplType = haplType
        bp.haplIdx = haplIdx
        bp.position = position
        bp.pairedPosition = pairedPosition
        bp.insertStr = insertStr

        return bp
