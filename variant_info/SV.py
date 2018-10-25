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
        self.copyNumber = -1
        self.genotype = ""  # P0P1M0
        self.haplRemain = None

    def info_str_title(self):
        return "length\tcopy_number\tgenotype\thapl_remain"

    def info_str(self):
        return "{0}\t{1}\t{2}\t{3}".format(
            self.length, self.copyNumber, self.genotype, self.haplRemain)


class INSERTION:

    def __init__(self):

        self.length = -1
        self.haplType = ""
        self.haplIdx = -1

    def info_str_title(self):
        return "length\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}".format(
            self.length, self.haplType, self.haplIdx)


class DELETION:

    def __init__(self):

        self.length = -1
        self.haplType = ""
        self.haplIdx = -1

    def info_str_title(self):
        return "length\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}".format(
            self.length, self.haplType, self.haplIdx)


class COMPLEXINDEL:

    def __init__(self):

        self.length1 = -1
        self.length2 = -1
        self.haplType = ""
        self.haplIdx = -1

    def info_str_title(self):
        return "length\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}".format(
            self.length, self.haplType, self.haplIdx)


class INVERTION:

    def __init__(self):

        self.length = -1
        self.haplType = ""
        self.haplIdx = -1

    def info_str_title(self):
        return "length\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}".format(
            self.length, self.haplType, self.haplIdx)


class TANDEMDUP:

    def __init__(self):

        self.length = -1
        self.haplType = ""
        self.haplIdx = -1
        self.times = 1

    def info_str_title(self):
        return "length\times\thapl_type\thapl_idx"

    def info_str(self):
        return "{0}\t{1}\t{2}\t{3}".format(
            self.length, self.times, self.haplType, self.haplIdx)


class TRANSLOCATION:

    def __init__(self):
        self.chromTo = ""
        self.haplTypeFrom = ""
        self.haplTypeTo = ""
        self.haplIdxFrom = -1
        self.haplIdxTo = -1
        self.length = -1

    def info_str_title(self):
        return "chromTo\thapl_type_from\thapl_type_to\
            \thapl_idx_from\thapl_idx_to\tlength"

    def info_str(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(
            self.chromTo, self.haplTypeFrom, self.haplTypeTo,
            self.haplIdxFrom, self.haplIdxTo, self.length)



class OVERLAPPEDCNV:


    def __init__(self):

        self.chrom = ""
        self.haplType = ""
        self.haplIdx = -1
        self.length = -1
        self.inPositionStart = -1
        self.inPositionEnd = -1
        self.copyAddedNumber = 0

    def info_str_title(self):
        return "chrom\thaplType\thaplIdx\tinPositionStart\tinPositionEnd\tlength\tcopy_added_number"

    def info_str(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}".format(self.chrom, self.haplType,
                                                self.haplIdx,
                                                self.inPositionStart,
                                                self.inPositionEnd, self.length,
                                                self.copyAddedNumber)


class RTRANSLOCATION:

    def __init__(self):
        self.chromFrom = ""
        self.haplTypeFrom = ""
        self.haplIdxFrom = -1
        self.positionFrom = -1
        self.chromTo = ""
        self.haplTypeTo = ""
        self.haplIdxTo = -1
        self.positionTo = -1

    def info_str_title(self):
        return "chrom_from\thapl_type_from\thapl_idx_from\
            \tchrom_to\thapl_type_to\hapl_idx_to"

    def info_str(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(
            self.chromFrom, self.haplTypeFrom, self.haplIdxFrom,
            self.chromTo, self.haplTypeTo, self.haplIdxTo)

