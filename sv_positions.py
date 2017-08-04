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
from variant_info.SV import CNV, INDEL, INVERSION, TRANSVERTION


class SVP(object):

    """Docstring for SVP. """

    def __init__(self):

        # [variant_length, variant_copy_number,
        # variant_genotype, variant_name])
        self.position = -1
        self.sv_type = ""  # CNV, INDEL, TRANSVERSION, INVERSION

        self.sv = None


class SV_positions:

    def __init__(self):

        self.svp_dict = {}

    def add_posi_CNV(self, chrom, position, length, copy_number, genotype):
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

        self._add_svp(chrom, temp_SVP)

    def add_posi_INDEL(self, chrom, position, hapl_key, length, indel_type):
        temp_SVP = SVP()
        temp_SVP.sv_type = "INDEL"
        temp_SVP.position = position

        temp_SVP.sv = INDEL()
        temp_SVP.sv.indel_type = indel_type
        temp_SVP.sv.hapl_key = hapl_key
        temp_SVP.sv.length = length

        self._add_svp(chrom, temp_SVP)

    def add_posi_INVERSION(self, chrom, position, hapl_key, length):
        temp_SVP = SVP()
        temp_SVP.sv_type = "INVERSION"
        temp_SVP.position = position

        temp_SVP.sv = INVERSION()
        temp_SVP.sv.hapl_key = hapl_key
        temp_SVP.sv.length = length

        self._add_svp(chrom, temp_SVP)

    def add_posi_TRANSVERSION(
            self,
            chrom,
            position,
            hapl_key_from,
            hapl_key_to,
            length):
        temp_SVP = SVP()
        temp_SVP.sv_type = "TRANSVERSION"
        temp_SVP.position = position

        temp_SVP.sv = TRANSVERSION()
        temp_SVP.sv.hapl_key_from = hapl_key_from
        temp_SVP.sv.hapl_key_to = hapl_key_to
        temp_SVP.sv.length = length

        self._add_svp(chrom, temp_SVP)

    # 此处似乎不需要
    # def add_ploidy(self, chrom, hapl, number):


    def _add_svp(self, chrom, svp):
        if chrom not in svp_dict.keys():
            self.svp_dict[chrom] = [svp]
        else:
            self.svp_dict[chrom].append(svp)

    def sorted(self):
        for chrom in self.svp_dict.keys():
            self.svp_dict[chrom] = sorted(
                self.svp_dict[chrom],
                key=lambda d: d.position, reverse=True)
