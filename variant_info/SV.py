#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: SV.py
#          Desc: SV info
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-08-03 10:22:34
#       History:
# =============================================================================
'''


class CNV:

    def __init__(self):

        self.length = -1
        self.copy_number = -1
        self.genotype = ""  # P0P1M0
        self.hapl_remain = None

    def info_str_title(self):
        return "length\tcopy_number\tgenotype\thapl_remain"

    def info_str(self):
        return "{0}\t{1}\t{2}\t{3}".format(
            self.length, self.copy_number, self.genotype, self.hapl_remain)


class INSERTION:

    def __init__(self):

        self.length = -1
        self.hapl_type = ""
        self.hapl_idx = -1

    def info_str_title(self):
        return "length\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}".format(
            self.length, self.hapl_type, self.hapl_idx)


class DELETION:

    def __init__(self):

        self.length = -1
        self.hapl_type = ""
        self.hapl_idx = -1

    def info_str_title(self):
        return "length\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}".format(
            self.length, self.hapl_type, self.hapl_idx)


class INVERTION:

    def __init__(self):

        self.length = -1
        self.hapl_type = ""
        self.hapl_idx = -1

    def info_str_title(self):
        return "length\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}".format(
            self.length, self.hapl_type, self.hapl_idx)


class TANDEMDUP:

    def __init__(self):

        self.length = -1
        self.hapl_type = ""
        self.hapl_idx = -1
        self.times = 1

    def info_str_title(self):
        return "length\times\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}\t{3}".format(
            self.length, self.times, self.hapl_type, self.hapl_idx)


class TRANSLOCATION:

    def __init__(self):
        self.chrom_to = ""
        self.hapl_type_from = ""
        self.hapl_type_to = ""
        self.hapl_idx_from = -1
        self.hapl_idx_to = -1
        self.length = -1

    def info_str_title(self):
        return "chrom_to\thapl_type_from\thapl_type_to\
            \thapl_idx_from\thapl_idx_to\tlength"

    def info_str(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(
            self.chrom_to, self.hapl_type_from, self.hapl_type_to,
            self.hapl_idx_from, self.hapl_idx_to, self.length)
