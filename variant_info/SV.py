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


class INDEL:

    def __init__(self):

        self.indel_type = ""
        self.hapl_key = ""
        self.length = -1


class INVERTION:

    def __init__(self):

        self.hapl_key = ""
        self.length = -1


class TRANSVERTION:

    def __init__(self):

        self.hapl_key_from = ""
        self.hapl_key_to = ""
        self.length = -1
