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

        self.B_allele = ""


class SNV_positions:

    def __init__(self):

        self.snvp_dict = {}

    def add_ploidy(self, chrom, hapl, hapl_ids):
        hic = Counter(hapl_ids)
        ploidies = []
        for hi in hic.keys():
            hapl_idex = int(hi)
            number = hic[hi]
            ploidies = ploidies +\
                number * [self.snvp_dict[chrom][hapl][hapl_idex]]

        self.snvp_dict[chrom][hapl] = self.snvp_dict[chrom][hapl] + ploidies

        # self.breakpoints.delete_ploidy(chrom, hapl, hapl_idxes)
    def delete_ploidy(self, chrom, hapl, hapl_idxes):
        self.snvp_dict[chrom][hapl] = [
            self.snvp_dict[chrom][hapl][i]
            for i in range(self.snvp_dict[chrom][hapl])
            if i not in hapl_idxes]


    def addPosi(
            self,
            chrom,
            hapl_type,
            hapl_index,
            position,
            B_allele,
            isHetero,
            isOverlap):

        snvp = SNVP()
        snvp.position = position
        snvp.isHetero = isHetero
        snvp.isOverlap = isOverlap
        snvp.B_allele = B_allele

        self._add_snvp(chrom, hapl_type, hapl_index, snvp)

    def _add_snvp(self, chrom, hapl_type, hapl_index, snvp):

        if chrom not in self.snvp_dict.keys():
            self.snvp_dict[chrom] = {}
        if hapl_type not in self.snvp_dict[hapl_type].keys():
            self.snvp_dict[chrom][hapl_type] = []
        if hapl_index >= len(self.snvp_dict[chrom][hapl_type]):
            for i in range(hapl_index - len(self.snvp_dict[chrom][hapl_type])):
                self.snvp_dict[chrom][hapl_type].append([])

        self.snvp_dict[chrom][hapl_type][hapl_index].append(snvp)

    def sorted(self):

        for chrom in self.snvp_dict.keys():
            for hapl_type in self.snvp_dict[chrom].keys():
                for hapl_index in range(len(self.snvp_dict[chrom][hapl_type])):
                    self.snvp_dict[chrom][hapl_type] = sorted(
                        self.snvp_dict[chrom][hapl_type],
                        key=lambda item: item.position, reverse=True)
