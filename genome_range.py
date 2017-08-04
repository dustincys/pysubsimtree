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

from utils import compute_range


class GenomeRange(object):

    """GenomeRange, record positions"""

    def __init__(self, ref):
        """initialized by ref"""
        self._genome_range = self._init_from_ref(ref)

    def _init_from_ref(self, ref):
        return compute_range(ref)

    def sample1posi(self, chrom):
        return random.sample(self._genome_range[chrom], 1)[0]

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
        self._genome_range[chrom] = np.delete(self._genome_range[chrom],
                                              temp_range)
