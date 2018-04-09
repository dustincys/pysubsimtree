#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: genome_range.py
#          Desc: gnome position operations
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-07-28 08:36:41
#       History:
# =============================================================================
'''

import numpy as np
import random

from utils.utils import compute_range


class GenomeRange(object):

    """GenomeRange, record positions"""

    def __init__(self):
        """initialized by ref"""
        self._genomeRange = {}

    def init_none(self, ref):
        for chrom in ref.keys():
            self._genomeRange[chrom] = np.arange(0, 0)

    def init_from_ref(self, ref):
        # {'chr1':[], 'chr2':[] }
        self._genomeRange = compute_range(ref)

    def sample1posi(self, chrom):
        return random.sample(self._genomeRange[chrom], 1)[0]

    def samplePosis(self, chrom, number):
        return random.sample(self._genomeRange[chrom], number)

    def isOverlaped(self, chrom, start, end):
        tempRange = np.arange(start, end)
        if np.array_equal(np.intersect1d(self._genomeRange[chrom], tempRange),
                          tempRange):
            return True
        else:
            return False

    def takePosi(self, chrom, start, end):
        if not self.isOverlaped(chrom, start, end):
            print "Warning! remove not overlapping position!"
        tempRange = np.arange(start, end)
        rmIndex = np.where(np.in1d(self._genomeRange[chrom], tempRange))[0]
        self._genomeRange[chrom] = np.delete(self._genomeRange[chrom],
                                              rmIndex)

    def addRange(self, chrom, start, end):
        preRange = set(self._genomeRange[chrom])
        addRange = set(np.arange(start, end))
        self._genomeRange[chrom] = np.array(list(preRange | addRange))
