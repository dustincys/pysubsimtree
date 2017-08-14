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
        self._genome_range = {}

    def init_none(self, ref):
        for chrom in ref.keys():
            self._genome_range[chrom] = np.arange(0, 0)

    def init_from_ref(self, ref):
        # {'chr1':[], 'chr2':[] }
        self._genome_range = compute_range(ref)

    def sample1posi(self, chrom):
        return random.sample(self._genome_range[chrom], 1)[0]

    def samplePosis(self, chrom, number):
        return random.sample(self._genome_range[chrom], number)

    def isOverlaped(self, chrom, start, end):
        temp_range = np.arange(start, end)
        if np.array_equal(np.intersect1d(self._genome_range[chrom], temp_range),
                          temp_range):
            return True
        else:
            return False

    def takePosi(self, chrom, start, end):
        if not self.isOverlaped(chrom, start, end):
            print "Warning! remove not overlapping position!"
        temp_range = np.arange(start, end)
        rm_index = np.where(np.in1d(self._genome_range[chrom], temp_range))[0]
        self._genome_range[chrom] = np.delete(self._genome_range[chrom],
                                              rm_index)

    def addRange(self, chrom, start, end):
        pre_range = set(self._genome_range[chrom])
        add_range = set(np.arange(start, end))
        self._genome_range[chrom] = np.array(list(pre_range | add_range))
