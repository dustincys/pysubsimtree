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
        self.sv_type = ""  # CNV, INDEL, TRANSLOCATION, INVERSION

        self.sv = None

    def info_str_title(self):
        return "# position\tsv_type\t{}\n".format(self.sv.info_str_title())

    def info_str(self):
        return "{0}\t{1}\t{2}\n".format(self.position, self.sv_type,
                                        self.sv.info_str())


class SV_positions:

    def __init__(self):

        self.svp_dict = {}

    def add_posi_CNV(self, chrom, position, length, copy_number, genotype,
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
        temp_SVP = SVP()
        temp_SVP.sv_type = "CNV"
        temp_SVP.position = position

        temp_SVP.sv = CNV()
        temp_SVP.sv.length = length
        temp_SVP.sv.copy_number = copy_number
        temp_SVP.sv.genotype = genotype

        temp_SVP.sv.hapl_remain = self._getHaplRemain(
            genotype, ploidy_status[chrom])

        self._add_svp(chrom, temp_SVP)

    def _getHaplRemain(self, genotype, ploidy_status):
        if genotype == "NONE":
            return {}

        hr = {}

        gtc = Counter(genotype)
        gpsc = Counter(ploidy_status)

        for hapl_type in set(gtc.keys()) & set(gpsc.keys()):
            number = 0
            if gtc[hapl_type] >= gpsc[hapl_type]:
                number = gpsc[hapl_type]
            else:
                number = gtc[hapl_type]

            hr[hapl_type] = random.sample(range(gpsc[hapl_type]), number)

        return hr

    def add_posi_INSERTION(self, chrom, hapl_type, hapl_idx, position, length):
        temp_SVP = SVP()
        temp_SVP.sv_type = "INSERTION"
        temp_SVP.position = position

        temp_SVP.sv = INSERTION()
        temp_SVP.sv.hapl_type = hapl_type
        temp_SVP.sv.hapl_idx = hapl_idx
        temp_SVP.sv.length = length

        self._add_svp(chrom, temp_SVP)

    def add_posi_DELETION(self, chrom, hapl_type, hapl_idx, position, length):
        temp_SVP = SVP()
        temp_SVP.sv_type = "DELETION"
        temp_SVP.position = position

        temp_SVP.sv = DELETION()
        temp_SVP.sv.hapl_type = hapl_type
        temp_SVP.sv.hapl_idx = hapl_idx
        temp_SVP.sv.length = length

        self._add_svp(chrom, temp_SVP)

    def add_posi_COMPLEXINDEL(self, chrom, hapl_type, hapl_idx, position, length1, length2):
        temp_SVP = SVP()
        temp_SVP.sv_type = "COMPLEXINDEL"
        temp_SVP.position = position

        temp_SVP.sv = COMPLEXINDEL()
        temp_SVP.sv.hapl_type = hapl_type
        temp_SVP.sv.hapl_idx = hapl_idx
        temp_SVP.sv.length1 = length1
        temp_SVP.sv.length2 = length2

        self._add_svp(chrom, temp_SVP)

    def add_posi_INVERTION(self, chrom, hapl_type, hapl_idx, position, length):
        temp_SVP = SVP()
        temp_SVP.sv_type = "INVERTION"
        temp_SVP.position = position

        temp_SVP.sv = INVERTION()
        temp_SVP.sv.hapl_type = hapl_type
        temp_SVP.sv.hapl_idx = hapl_idx
        temp_SVP.sv.length = length

        self._add_svp(chrom, temp_SVP)

    def add_posi_TANDEMDUP(self, chrom, hapl_type, hapl_idx, position, length,
                           times):
        temp_SVP = SVP()
        temp_SVP.sv_type = "TANDEMDUP"
        temp_SVP.position = position

        temp_SVP.sv = TANDEMDUP()
        temp_SVP.sv.hapl_type = hapl_type
        temp_SVP.sv.hapl_idx = hapl_idx
        temp_SVP.sv.length = length
        temp_SVP.sv.times = times

        self._add_svp(chrom, temp_SVP)

    def add_posi_TRANSLOCATION(
            self,
            chrom_from,
            position_from,
            hapl_type_from,
            hapl_idx_from,
            chrom_to,
            hapl_type_to,
            hapl_idx_to,
            length):

        temp_SVP = SVP()
        temp_SVP.sv_type = "TRANSLOCATION"
        temp_SVP.position = position_from

        temp_SVP.sv = TRANSLOCATION()

        temp_SVP.sv.chrom_to = chrom_to

        temp_SVP.sv.hapl_type_from = hapl_type_from
        temp_SVP.sv.hapl_type_to = hapl_type_to
        temp_SVP.sv.hapl_idx_from = hapl_idx_from
        temp_SVP.sv.hapl_idx_to = hapl_idx_to
        temp_SVP.sv.length = length

        self._add_svp(chrom_from, temp_SVP)

    # 此处似乎不需要
    # def add_ploidy(self, chrom, hapl, number):

    def _add_svp(self, chrom, svp):
        if chrom not in self.svp_dict.keys():
            self.svp_dict[chrom] = [svp]
        else:
            self.svp_dict[chrom].append(svp)

    def sorted(self):
        for chrom in self.svp_dict.keys():
            self.svp_dict[chrom] = sorted(
                self.svp_dict[chrom],
                key=lambda d: d.position, reverse=True)
