#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: snv_positions.py
#          Desc: snv
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-08-07 14:54:48
#       History:
# =============================================================================
'''
from collections import Counter


class SNVP(object):

    """Docstring for SVP. """

    def __init__(self):

        self.position = -1

        self.isHetero = ""
        self.isOverlap = ""

        self.BAllele = ""


class SNVPositions:

    def __init__(self):

        self.snvpDict = {}

    def add_ploidy(self, chrom, hapl, haplIds):
        hic = Counter(haplIds)
        ploidies = []
        for hi in hic.keys():
            haplIdx = int(hi)
            number = hic[hi]
            if self._has(chrom, hapl, haplIdx):
                ploidies = ploidies +\
                    number * [self.snvpDict[chrom][hapl][haplIdx]]

        if self._has(chrom, hapl, haplIdx):
            self.snvpDict[chrom][hapl] = self.snvpDict[chrom][hapl] + ploidies

    def _has(self, chrom, hapl, haplIdx):
        if chrom in self.snvpDict.keys():
            if hapl in self.snvpDict[chrom].keys():
                if haplIdx < len(self.snvpDict[chrom][hapl]):
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def delete_ploidy(self, chrom, hapl, haplIdxes):
        self.snvpDict[chrom][hapl] = [
            self.snvpDict[chrom][hapl][i]
            for i in range(len(self.snvpDict[chrom][hapl]))
                if i not in haplIdxes]

    def addPosi(
            self,
            chrom,
            haplType,
            haplIndex,
            position,
            BAllele,
            isHetero,
            isOverlap):

        snvp = SNVP()
        snvp.position = position
        snvp.isHetero = isHetero
        snvp.isOverlap = isOverlap
        snvp.BAllele = BAllele

        self._add_snvp(chrom, haplType, haplIndex, snvp)

    def _add_snvp(self, chrom, haplType, haplIndex, snvp):

        if chrom not in self.snvpDict.keys():
            self.snvpDict[chrom] = {}
        if haplType not in self.snvpDict[chrom].keys():
            self.snvpDict[chrom][haplType] = []
        if haplIndex >= len(self.snvpDict[chrom][haplType]):
            for i in range(haplIndex+1-len(self.snvpDict[chrom][haplType])):
                self.snvpDict[chrom][haplType].append([])

        self.snvpDict[chrom][haplType][haplIndex].append(snvp)

    def sorted(self):

        for chrom in self.snvpDict.keys():
            for haplType in self.snvpDict[chrom].keys():
                for haplIndex in range(len(self.snvpDict[chrom][haplType])):
                    self.snvpDict[chrom][haplType] = sorted(
                        self.snvpDict[chrom][haplType],
                        key=lambda item: item.position, reverse=True)
