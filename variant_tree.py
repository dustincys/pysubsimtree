#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: genomeNodes.py
#          Desc: 根据配置文件，生成变异树，根据要求可以遍历树的每一个节点
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-07-25 22:45:00
#       History:
# =============================================================================
'''

from breakpoint import BreakPoints
from sv_positions import SVPositions
from snv_positions import SNVPositions
from genome_range import GenomeRange

from utils.utils import outputFa, SNP
from utils import constants

from anytree import NodeMixin

import copy
from collections import Counter
import random


genomeNodeType = ["SV", "SNV", "PLOIDY"]


class GenomeNode(NodeMixin):
    def __init__(self, name, parent=None, variantType="SV"):
        assert variantType in genomeNodeType

        self.name = name
        self.parent = parent
        # 用来传递生成ＳＮＶ

        self.variantType = variantType

        self.variantL = []

        # generate ploidy reference
        self.ploidyStatus = {}

        self.snvPositions = None
        self.svPositions = None
        self.breakpoints = None
        self.availPosition = None

        # record the none del sv position
        self.noDELPosition = None

    def outputRef(self, ref, outfilePrefix):
        refPloidy = self._ref_make_ploidy(copy.deepcopy(ref))
        refPloidySNV = self._ref_make_snv(refPloidy)
        variantRef = self.breakpoints.generateFa(refPloidySNV)
        outputFa(variantRef, outfilePrefix+self.name+".fa")

    def _ref_make_ploidy(self, ref):
        refPloidy = {}
        for chrom in self.ploidyStatus.keys():
            if chrom not in refPloidy.keys():
                refPloidy[chrom] = {}

            psc = Counter(self.ploidyStatus[chrom])
            for haplType in psc.keys():
                if haplType not in refPloidy[chrom].keys():
                    refPloidy[chrom][haplType] = []

                for i in range(psc[haplType]):
                    refPloidy[chrom][haplType].append(
                        ref[chrom][haplType][0])

        return refPloidy

    def _ref_make_snv(self, refPloidy):
        # 与breakpoint结构一致
        for chrom in self.snvPositions.snvpDict.keys():
            for haplType in self.snvPositions.snvpDict[chrom].keys():
                for haplIdx in \
                        range(len(self.snvPositions.snvpDict
                                  [chrom][haplType])):
                    lstr = list(refPloidy[chrom][haplType][haplIdx])
                    for snvp in\
                            self.snvPositions.snvpDict[
                                                chrom][haplType][haplIdx]:
                        lstr[snvp.position] = snvp.BAllele

                    refPloidy[chrom][haplType][haplIdx] = ''.join(lstr)

        return refPloidy

    def outputVartPosi(self, outfilePrefix):
        # self.svPositions[chroms[k]].append([
        # posi,
        # self.variantL[i][0],
        # self.variantL[i][1],
        # self.variantL[i][2],
        # self.variantL[i][-1]])
        # tempNode.variantL.append(
        # [variantLength, variantCopyNumber,
        # variantGenotype, variantName])
        outfileName = outfilePrefix+self.name+".SV.txt"

        with open(outfileName, 'w') as outfile:
            for chrom in self.svPositions.svpDict.keys():
                svItems = self.svPositions.svpDict[chrom]
                for svItem in svItems:
                    outfile.write("chrom\t" + svItem.info_str_title())
                    outfile.write(chrom + "\t" + svItem.info_str())
        pass

    def make(self, ref):
        """generate variats
        :returns: TODO

        """
        # 首先初始化节点, 即从父节点继承信息
        self._init_node(ref)

        # 第一步，生成ploidy 变异
        if self.variantType == "PLOIDY":
            self._make_ploidy()
        else:

            # 第一步，生成非结构变异，是否与结构变异出现重合现象？
            # 生成的位置包括： available中位置， insertion中的位置

            refPloidy = self._ref_make_ploidy(copy.deepcopy(ref))
            refPloidySNV = self._ref_make_snv(refPloidy)

            # 第二步，生成结构变异, 此处是生成了ploidy的ref
            self._make_sv(refPloidySNV)
            self._make_snv(refPloidySNV)

    def _init_node(self, ref):
        self.availPosition = self._getPrtAvailPst(ref)
        self.svPositions = self._getPrtSVPosis()
        self.snvPositions = self._getPrtSNVPosis()
        self.breakpoints = self._getPrtBreakpoints()
        self.ploidyStatus = self._getPrtPloidyStatus(ref)
        self.noDELPosition = self._getPrtNBPPst(ref)

    def _make_ploidy(self):
        # 此处对从父节点继承而来的变异信息进行染色体复制操作
        # 单体型需要用字典来实现ploidy变异
        # 从父节点的breakpoints和svPositions中生成

        # variant = [chrom, ploidyNumberBefore,
                    # ploidyNumberAfter, ploidyTypeBefore,
                    # ploidyTypeAfter]

        for variant in self.variantL:
            chrom = variant[0]
            # ploidyNumberBefore = variant[1]
            # ploidyNumberAfter = variant[2]
            ploidyTypeBefore = variant[3]
            ploidyTypeAfter = variant[4]
            ptcBefore = Counter(ploidyTypeBefore)
            ptcAfter = Counter(ploidyTypeAfter)

            for hapl in ptcBefore.keys():
                if hapl in ptcAfter.keys():
                    if ptcBefore[hapl] > ptcAfter[hapl]:
                        # 此处可能不需要针对sv_position进行ploidy
                        # 操作，而只需要对breakpoint 进行操作
                        number = ptcBefore[hapl]-ptcAfter[hapl]
                        haplIdxes = random.sample(
                            range(ptcBefore[hapl]), number)
                        self.breakpoints.delete_ploidy(chrom, hapl, haplIdxes)
                        self.snvPositions.delete_ploidy(chrom, hapl, haplIdxes)

                    elif ptcBefore[hapl] < ptcAfter[hapl]:
                        number = ptcAfter[hapl]-ptcBefore[hapl]
                        hapls = [random.choice(range(ptcBefore[hapl])) for _ in
                                 range(number)]
                        self.breakpoints.add_ploidy(chrom, hapl, hapls)
                        self.snvPositions.add_ploidy(chrom, hapl, hapls)

                    else:
                        continue
                else:
                    number = ptcBefore[hapl]
                    haplIdxes = range(ptcBefore[hapl])
                    self.breakpoints.delete_ploidy(chrom, hapl, haplIdxes)
                    self.snvPositions.delete_ploidy(chrom, hapl, haplIdxes)

        self.ploidyStatus[chrom] = ploidyTypeAfter
        pass

    def _make_snv(self, ref):
        # 从availPosition 中生成非重合变异，
        # 从sv_position和breakpoints的insertStr中生成重合的变异
        # 只有非重合变异可以生成纯合突变

        snvList = filter(lambda item: item[-1] == "SNV", self.variantL)
        for snv in snvList:
            chrom = snv[0]
            isHetero = snv[1]
            isOverlap = snv[2]
            number = snv[3]

            if isOverlap == "TRUE":
                # only hetero
                numberNoDEL = number
                numberBp = 0
                noDELCNVpois = self._getNDCPs(chrom)
                strBps = self._getStrBps(chrom)

                if strBps is not None and len(noDELCNVpois) > 0:
                    numberNoDEL = int(number / 2)
                    numberBp = number - numberNoDEL
                if len(noDELCNVpois) == 0:
                    numberNoDEL = 0
                    numberBp = number
                if strBps is None:
                    numberBp = 0
                    numberNoDEL = number
                if strBps is None and len(noDELCNVpois) == 0:
                    numberBp = 0
                    numberNoDEL = 0

                for i in range(numberBp):
                    bp = random.sample(strBps, 1)[0]
                    position = random.sample(range(len(bp.insertStr)), 1)[0]
                    bpLstr = list(bp.insertStr)
                    bpLstr[position] = SNP(bpLstr[position])
                    bp.insertStr = ''.join(bpLstr)

                for i in range(numberNoDEL):

                    noDELCNV = random.sample(noDELCNVpois, 1)[0]
                    length = noDELCNV.sv.length
                    haplRemain = noDELCNV.sv.haplRemain
                    position = noDELCNV.position\
                        + random.sample(range(length), 1)[0]

                    tempHaplType = random.sample(haplRemain.keys(), 1)[0]
                    tempHaplIndex = random.sample(
                        haplRemain[tempHaplType], 1)[0]
                    haplLstr = ref[chrom][tempHaplType][tempHaplIndex]
                    BAllele = SNP(haplLstr[position])

                    self.snvPositions.addPosi(chrom, tempHaplType,
                                               tempHaplIndex,
                                               position, BAllele, isHetero,
                                               isOverlap)
            else:

                positions = self.availPosition.samplePosis(chrom, number)

                if isHetero == "TRUE":
                    haplType = random.sample(ref[chrom].keys(), 1)[0]

                    for pi in positions:
                        haplIndex = random.sample(
                            range(len(ref[chrom][haplType])), 1)[0]
                        haplLstr = ref[chrom][haplType][haplIndex]
                        BAllele = SNP(haplLstr[pi])

                        self.snvPositions.addPosi(chrom, haplType, haplIndex,
                                                   pi, BAllele, isHetero,
                                                   isOverlap)
                else:
                    for pi in positions:
                        tempHaplType = random.sample(ref[chrom].keys(), 1)[0]
                        tempHaplIndex = random.sample(
                            range(len(ref[chrom][tempHaplType])), 1)[0]
                        BAllele = SNP(
                            ref[chrom][tempHaplType][tempHaplIndex][pi]
                        )
                        for haplType in ref[chrom].keys():
                            for haplIndex in range(len(ref[chrom][haplType])):
                                self.snvPositions.addPosi(chrom, haplType,
                                                           haplIndex,
                                                           pi, BAllele,
                                                           isHetero, isOverlap)
        pass

    def _getStrBps(self, chrom):
        # 此处为了进行ploidy操作，bp结构为 {chr1:{'P':[[bp1],[bp2]], 'M':[[bp3]]}}

        for haplType in self.breakpoints.breaksList[chrom].keys():
            for haplIndex in range(len(
                    self.breakpoints.breaksList[chrom][haplType])):
                strBps = filter(lambda item: item.insertStr is not None,
                                self.breakpoints.breaksList
                                [chrom][haplType][haplIndex])
                if len(strBps) > 1:
                    return strBps

        return None

    def _getNDCPs(self, chrom):
        return filter(
            lambda item: item.svType == "CNV" and item.sv.genotype != "NONE",
            self.svPositions.svpDict[chrom])

    def _make_sv(self, ref):
        newPositions = self._generateNewPosition()
        # noBP remove deletion seg

        # 此处对只断点内部所复制的字符串进行操作，且断点内部的字符串是由父节点生
        # 成
        # 此处应该生成新位置，但不生成新断点
        self.breakpoints.generateOverlappedCNV(newPositions)

        self.breakpoints.generateBPs(newPositions, self.availPosition,
                                     self.noDELPosition, ref,
                                     self.ploidyStatus)


    def _generateNewPosition(self):
        clmin = constants.COMPLEXINDEL_LENGTH_MIN
        clmax = constants.COMPLEXINDEL_LENGTH_MAX

        currentSVPositions = SVPositions()

        # str_len=[len(ref[tmp][0]) for tmp in dic]
        svList = filter(lambda item: item[-1] == "SV", self.variantL)

        for sv in svList:
            # 此处sv变异被按顺序分配到染色体上
            chrom = sv[0]
            length = sv[1]
            variantName = sv[-2]


            # 此处添加重叠CNV变异
            if variantName == "OVERLAPPEDCNV":
                position = sv[2]
                inPositionStart = sv[3]
                inPositionEnd = sv[4]
                haplType = sv[5]
                haplIdx = sv[6]
                copyAddedNumber = sv[7]

                self.svPositions.add_posi_OVERLAPPEDCNV(
                    chrom, position, inPositionStart, inPositionEnd, haplType, haplIdx, length, copyAddedNumber)
                currentSVPositions.add_posi_OVERLAPPEDCNV(
                    chrom, position, inPositionStart, inPositionEnd, haplType, haplIdx, length, copyAddedNumber)
                continue

            while True:
                count = 1
                posi = self.availPosition.sample1posi(chrom)
                # 生成位置，以随机生成的第一个首字母为起始位置。
                if variantName == "INSERTION":
                    isOverlap = self.availPosition.isOverlaped(chrom, posi,
                                                                posi + 1)
                else:
                    isOverlap = self.availPosition.isOverlaped(chrom, posi,
                                                                posi + length)

                if isOverlap:
                    # 此处需要在conf 中给定chrom, 因为不同的染色体的ploidy
                    # 不同，对应的变异的基因型也不同
                    if variantName == "INSERTION":
                        self.availPosition.takePosi(chrom, posi, posi + 1)
                    else:
                        self.availPosition.takePosi(chrom, posi, posi + length)

                    if variantName == "CNV":
                        copyNumber = sv[2]
                        genotype = sv[3]

                        self.svPositions.add_posi_CNV(
                            chrom, posi, length, copyNumber, genotype,
                            self.ploidyStatus)
                        currentSVPositions.add_posi_CNV(
                            chrom, posi, length, copyNumber, genotype,
                            self.ploidyStatus)

                    if variantName == "INVERTION":

                        haplType = sv[2]
                        haplIdx = sv[3]
                        # ploidyStatus = sv[4]

                        self.svPositions.add_posi_INVERTION(
                            chrom, haplType, haplIdx, posi, length)
                        currentSVPositions.add_posi_INVERTION(
                            chrom, haplType, haplIdx, posi, length)

                    if variantName == "TANDEMDUP":

                        haplType = sv[2]
                        haplIdx = sv[3]
                        times = sv[4]
                        # ploidyStatus = sv[5]

                        self.svPositions.add_posi_TANDEMDUP(
                            chrom, haplType, haplIdx, posi, length, times)
                        currentSVPositions.add_posi_TANDEMDUP(
                            chrom, haplType, haplIdx, posi, length, times)

                    if variantName == "DELETION":

                        haplType = sv[2]
                        haplIdx = sv[3]
                        # ploidyStatus = sv[4]

                        self.svPositions.add_posi_DELETION(
                            chrom, haplType, haplIdx, posi, length)
                        currentSVPositions.add_posi_DELETION(
                            chrom, haplType, haplIdx, posi, length)

                    if variantName == "COMPLEXINDEL":

                        haplType = sv[2]
                        haplIdx = sv[3]
                        # ploidyStatus = sv[4]
                        length1 = random.randint(clmin, clmax)
                        length2 = random.randint(clmin, clmax)

                        self.svPositions.add_posi_COMPLEXINDEL(
                            chrom, haplType, haplIdx, posi, length1, length2)
                        currentSVPositions.add_posi_COMPLEXINDEL(
                            chrom, haplType, haplIdx, posi, length1, length2)

                    if variantName == "INSERTION":

                        haplType = sv[2]
                        haplIdx = sv[3]
                        # ploidyStatus = sv[4]

                        self.svPositions.add_posi_INSERTION(
                            chrom, haplType, haplIdx, posi, length)
                        currentSVPositions.add_posi_INSERTION(
                            chrom, haplType, haplIdx, posi, length)

                    if variantName == "TRANSLOCATION":
                        haplTypeFrom = sv[2]
                        haplIdxFrom = int(sv[3])
                        chromTo = sv[4]
                        haplTypeTo = sv[5]
                        haplIdxTo = int(sv[6])
                        # ploidy_genotype = sv[7]

                        self.svPositions.add_posi_TRANSLOCATION(
                            chrom, posi, haplTypeFrom, haplIdxFrom,
                            chromTo, haplTypeTo, haplIdxTo, length)

                        currentSVPositions.add_posi_TRANSLOCATION(
                            chrom, posi, haplTypeFrom, haplIdxFrom,
                            chromTo, haplTypeTo, haplIdxTo, length)
                    break
                else:
                    count = count+1
                    if count < 200:
                        continue
                    else:
                        break

        self.svPositions.sorted()
        currentSVPositions.sorted()

        return currentSVPositions

    def _getPrtBreakpoints(self):
        if self._isRoot():
            return BreakPoints()
        else:
            return copy.deepcopy(self.parent.breakpoints)
        pass

    def _getPrtSNVPosis(self):
        if self._isRoot():
            return SNVPositions()
        else:
            return copy.deepcopy(self.parent.snvPositions)
        pass

    def _getPrtSVPosis(self):
        if self._isRoot():
            return SVPositions()
        else:
            return copy.deepcopy(self.parent.svPositions)
        pass

    def _getPrtAvailPst(self, ref):
        if self._isRoot():
            availPosition = GenomeRange()
            availPosition.init_from_ref(ref)
            return availPosition
        else:
            return copy.deepcopy(self.parent.availPosition)

    def _getPrtNBPPst(self, ref):
        if self._isRoot():
            availPosition = GenomeRange()
            availPosition.init_none(ref)
            return availPosition
        else:
            return copy.deepcopy(self.parent.noDELPosition)

    def _getPrtPloidyStatus(self, ref):
        if self._isRoot():
            ps = {}
            for chrom in ref.keys():
                ps[chrom] = "PM"
            return ps
        else:
            return copy.deepcopy(self.parent.ploidyStatus)

    def _isRoot(self):
        if self.parent is None:
            return True
        else:
            return False


class GenomeTree(object):

    """变异树"""

    def __init__(self, configFileName):
        self.genomeNodes = self._read_config(configFileName)

    def _getParentName(self, name):
        """生成父节点的名称

        格式：
        1
        1.1
        1.2
        1.1.1
        ...

        :name: 当前节点名称
        :returns: 父节点名称

        """
        nameList = name.rsplit(".", 1)
        if len(nameList) == 1:
            parentName = None
        else:
            parentName = nameList[0]

        return parentName

    def _read_config(self, fileName):
        genomeNodes = {}

        with open(fileName) as infile:
            for line in sorted(infile):
                # 这里按顺序生成
                if not line.startswith('#'):
                    listLine = line.rstrip().split('\t')
                else:
                    continue

                subclonalName = listLine[0]
                variantType = listLine[1]

                if variantType == "SV":
                    variantName = listLine[2]

                    if variantName == "CNV":
                        chrom = listLine[3]
                        variantLength = int(listLine[4])
                        variantCopyNumber = int(listLine[5])
                        variantGenotype = listLine[6]
                        number = int(listLine[7])
                        variant = [
                            chrom,
                            variantLength,
                            variantCopyNumber,
                            variantGenotype,
                            variantName,
                            variantType]
                        self._add2node(
                            variant, number, subclonalName, genomeNodes)

                    elif variantName == "INVERTION":
                        chrom = listLine[3]
                        haplType = listLine[4]
                        haplIdx = int(listLine[5])
                        variantLength = int(listLine[6])
                        number = int(listLine[7])

                        variant = [
                            chrom,
                            variantLength,
                            haplType,
                            haplIdx,
                            variantName,
                            variantType]

                        self._add2node(
                            variant, number, subclonalName, genomeNodes)

                    elif variantName == "TANDEMDUP":
                        chrom = listLine[3]
                        haplType = listLine[4]
                        haplIdx = int(listLine[5])
                        variantLength = int(listLine[6])
                        times = int(listLine[7])
                        number = int(listLine[8])

                        variant = [
                            chrom,
                            variantLength,
                            haplType,
                            haplIdx,
                            times,
                            variantName,
                            variantType]

                        self._add2node(
                            variant, number, subclonalName, genomeNodes)

                    elif variantName == "INSERTION":
                        chrom = listLine[3]
                        haplType = listLine[4]
                        haplIdx = int(listLine[5])
                        variantLength = int(listLine[6])
                        number = int(listLine[7])

                        variant = [
                            chrom,
                            variantLength,
                            haplType,
                            haplIdx,
                            variantName,
                            variantType]

                        self._add2node(
                            variant, number, subclonalName, genomeNodes)

                    elif variantName == "DELETION":
                        chrom = listLine[3]
                        haplType = listLine[4]
                        haplIdx = int(listLine[5])
                        variantLength = int(listLine[6])
                        number = int(listLine[7])

                        variant = [
                            chrom,
                            variantLength,
                            haplType,
                            haplIdx,
                            variantName,
                            variantType]

                        self._add2node(
                            variant, number, subclonalName, genomeNodes)

                    elif variantName == "COMPLEXINDEL":
                        chrom = listLine[3]
                        haplType = listLine[4]
                        haplIdx = int(listLine[5])
                        number = int(listLine[6])

                        variant = [
                            chrom,
                            -1,
                            haplType,
                            haplIdx,
                            variantName,
                            variantType]

                        self._add2node(
                            variant, number, subclonalName, genomeNodes)

                    elif variantName == "TRANSLOCATION":
                        chrom = listLine[3]
                        haplTypeFrom = listLine[4]
                        haplIdxFrom = int(listLine[5])
                        variantLength = int(listLine[6])
                        chromTo = listLine[7]
                        haplTypeTo = listLine[8]
                        haplIdxTo = int(listLine[9])
                        number = int(listLine[10])

                        variant = [
                            chrom,
                            variantLength,
                            haplTypeFrom,
                            haplIdxFrom,
                            chromTo,
                            haplTypeTo,
                            haplIdxTo,
                            variantName,
                            variantType]

                        self._add2node(
                            variant, number, subclonalName, genomeNodes)

                    elif variantName == "OVERLAPPEDCNV":
                        # chrom, position, haplType, haplIdx, length, copyAddedNumber
                        chrom = listLine[3]
                        position = int(listLine[4])
                        inPositionStart = int(listLine[4])
                        inPositionEnd = int(listLine[5])
                        haplType = listLine[6]
                        haplIdx = int(listLine[7])
                        length = int(listLine[8])
                        copyAddedNumber = int(listLine[9])

                        variant = [
                            chrom,
                            length,
                            position,
                            inPositionStart,
                            inPositionEnd,
                            haplType,
                            haplIdx,
                            copyAddedNumber,
                            variantName,
                            variantType]

                        self._add2node(
                            variant, number, subclonalName, genomeNodes)

                    else:
                        pass

                elif variantType == "PLOIDY":
                    chrom = listLine[2]
                    ploidyNumberBefore = int(listLine[3])
                    ploidyNumberAfter = int(listLine[4])
                    ploidyTypeBefore = listLine[5]
                    ploidyTypeAfter = listLine[6]
                    variant = [chrom, ploidyNumberBefore,
                               ploidyNumberAfter, ploidyTypeBefore,
                               ploidyTypeAfter, variantType]

                    self._add2node(
                        variant, 1, subclonalName, genomeNodes)

                elif variantType == "SNV":
                    chrom = listLine[2]
                    isHetero = listLine[3]
                    isOverlap = listLine[4]
                    number = int(listLine[6])
                    variant = [chrom, isHetero, isOverlap, number,
                               variantType]
                    self._add2node(
                        variant, 1, subclonalName, genomeNodes)

        self._linkNode(genomeNodes)

        return genomeNodes

    def _add2node(self, variant, number, subclonalName, genomeNodes):
        if subclonalName not in genomeNodes.keys():
            tempNode = GenomeNode(subclonalName,
                                    variantType=variant[-1])
            genomeNodes[subclonalName] = tempNode
        else:
            tempNode = genomeNodes[subclonalName]

        for i in range(number):
            tempNode.variantL = tempNode.variantL + [variant]

    def _linkNode(self, genomeNodes):
        if "1" not in genomeNodes.keys():
            genomeNodes["1"] = GenomeNode("1")

        for nodeName in genomeNodes.keys():
            nodeParentName = self._getParentName(nodeName)
            if nodeParentName is not None:
                try:
                    genomeNodes[nodeName].parent =\
                        genomeNodes[nodeParentName]
                except KeyError:
                    print "Node parent error!"

    def makeSV(self, ref):
        self._make(self.genomeNodes['1'], ref)

    def _make(self, node, ref):
        node.make(ref)
        for child in node.descendants:
            self._make(child, ref)


def main():
    """
    :returns: TODO

    """
    # genomeNodes = read_config("../config/sv_config_test")
    # print RenderTree(genomeNodes['1'])
    # for pre, _, node in RenderTree(genomeNodes['1']):
    #     print node.name, node.variantL


if __name__ == "__main__":
    main()
