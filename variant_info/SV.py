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


class INSERTION:

    def __init__(self):

        self.length = -1
        self.hapl_type = ""
        self.hapl_idx = -1


class DELETION:

    def __init__(self):

        self.length = -1
        self.hapl_type = ""
        self.hapl_idx = -1


class INVERTION:

    def __init__(self):

        self.length = -1
        self.hapl_type = ""
        self.hapl_idx = -1


class TRANSLOCATION:

    def __init__(self):
        self.chrom_to = ""
        self.hapl_type_from = ""
        self.hapl_type_to = ""
        self.hapl_idx_from = -1
        self.hapl_idx_to = -1
        self.length = -1
