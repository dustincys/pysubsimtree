#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: sv_positions.py
#          Desc: sv
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-08-03 11:10:33
#       History:
# =============================================================================
'''
from variant_info.SV import CNV, INSERTION, DELETION, INVERTION, TRANSLOCATION,\
    TANDEMDUP, COMPLEXINDEL
from collections import Counter
import random


class SVP(object):

    """Docstring for SVP. """

    def __init__(self):

        # [variant_length, variant_copy_number,
        # variant_genotype, variant_name])
        self.position = -1
        self.svType = ""  # CNV, INDEL, TRANSLOCATION, INVERSION

        self.sv = None

    def info_str_title(self):
        return "# position\tsv_type\t{}\n".format(self.sv.info_str_title())

    def info_str(self):
        return "{0}\t{1}\t{2}\n".format(self.position, self.svType,
                                        self.sv.info_str())


class SVPositions:

    def __init__(self):

        self.svpDict = {}

    def add_posi_CNV(self, chrom, position, length, copyNumber, genotype,
                     ploidy_status):
        # if chroms[k] not in self.sv_positions.keys():
        # self.sv_positions[chroms[k]] = []
        # self.sv_positions[chroms[k]].append([
        # pois,
        # self.sv_list[i][0],
        # self.sv_list[i][1],
        # self.sv_list[i][2],
        # self.sv_list[i][-1]])
        # temp_Node.sv_list.append(
        # [variant_length, variant_copy_number,
        # variant_genotype, variant_name])
        tempSVP = SVP()
        tempSVP.svType = "CNV"
        tempSVP.position = position

        tempSVP.sv = CNV()
        tempSVP.sv.length = length
        tempSVP.sv.copyNumber = copyNumber
        tempSVP.sv.genotype = genotype

        tempSVP.sv.haplRemain = self._getHaplRemain(
            genotype, ploidy_status[chrom])

        self._add_svp(chrom, tempSVP)

    def _getHaplRemain(self, genotype, ploidy_status):
        if genotype == "NONE":
            return {}

        hr = {}

        gtc = Counter(genotype)
        gpsc = Counter(ploidy_status)

        for haplType in set(gtc.keys()) & set(gpsc.keys()):
            number = 0
            if gtc[haplType] >= gpsc[haplType]:
                number = gpsc[haplType]
            else:
                number = gtc[haplType]

            hr[haplType] = random.sample(range(gpsc[haplType]), number)

        return hr

    def add_posi_INSERTION(self, chrom, haplType, haplIdx, position, length):
        tempSVP = SVP()
        tempSVP.svType = "INSERTION"
        tempSVP.position = position

        tempSVP.sv = INSERTION()
        tempSVP.sv.haplType = haplType
        tempSVP.sv.haplIdx = haplIdx
        tempSVP.sv.length = length

        self._add_svp(chrom, tempSVP)

    def add_posi_DELETION(self, chrom, haplType, haplIdx, position, length):
        tempSVP = SVP()
        tempSVP.svType = "DELETION"
        tempSVP.position = position

        tempSVP.sv = DELETION()
        tempSVP.sv.haplType = haplType
        tempSVP.sv.haplIdx = haplIdx
        tempSVP.sv.length = length

        self._add_svp(chrom, tempSVP)

    def add_posi_COMPLEXINDEL(self, chrom, haplType, haplIdx, position, length1, length2):
        tempSVP = SVP()
        tempSVP.svType = "COMPLEXINDEL"
        tempSVP.position = position

        tempSVP.sv = COMPLEXINDEL()
        tempSVP.sv.haplType = haplType
        tempSVP.sv.haplIdx = haplIdx
        tempSVP.sv.length1 = length1
        tempSVP.sv.length2 = length2

        self._add_svp(chrom, tempSVP)

    def add_posi_INVERTION(self, chrom, haplType, haplIdx, position, length):
        tempSVP = SVP()
        tempSVP.svType = "INVERTION"
        tempSVP.position = position

        tempSVP.sv = INVERTION()
        tempSVP.sv.haplType = haplType
        tempSVP.sv.haplIdx = haplIdx
        tempSVP.sv.length = length

        self._add_svp(chrom, tempSVP)

    def add_posi_TANDEMDUP(self, chrom, haplType, haplIdx, position, length,
                           times):
        tempSVP = SVP()
        tempSVP.svType = "TANDEMDUP"
        tempSVP.position = position

        tempSVP.sv = TANDEMDUP()
        tempSVP.sv.haplType = haplType
        tempSVP.sv.haplIdx = haplIdx
        tempSVP.sv.length = length
        tempSVP.sv.times = times

        self._add_svp(chrom, tempSVP)

    def add_posi_TRANSLOCATION(
            self,
            chrF_from,
            position_from,
            haplTypeFrom,
            haplIdxFrom,
            chromTo,
            haplTypeTo,
            haplIdxTo,
            length):

        tempSVP = SVP()
        tempSVP.svType = "TRANSLOCATION"
        tempSVP.position = position_from

        tempSVP.sv = TRANSLOCATION()

        tempSVP.sv.chromTo = chromTo

        tempSVP.sv.haplTypeFrom = haplTypeFrom
        tempSVP.sv.haplTypeTo = haplTypeTo
        tempSVP.sv.haplIdxFrom = haplIdxFrom
        tempSVP.sv.haplIdxTo = haplIdxTo
        tempSVP.sv.length = length

        self._add_svp(chromFrom, tempSVP)

    # 此处似乎不需要
    # def add_ploidy(self, chrom, hapl, number):

    def _add_svp(self, chrom, svp):
        if chrom not in self.svpDict.keys():
            self.svpDict[chrom] = [svp]
        else:
            self.svpDict[chrom].append(svp)

    def sorted(self):
        for chrom in self.svpDict.keys():
            self.svpDict[chrom] = sorted(
                self.svpDict[chrom],
                key=lambda d: d.position, reverse=True)
