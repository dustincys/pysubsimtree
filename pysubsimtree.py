#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import time
import ConfigParser
import copy

from variant_tree import VariantTree
from utils.utils import read_pf_config, reference, read_dbsnp, generate_normal


def main():
        usage = """%prog -c ini.config Author: Yanshuo Chu"""
        parser = OptionParser(usage)
        parser.add_option(
            "-c",
            "--config",
            dest="config",
            help="config file",
            metavar="FILE")
        parser.add_option(
            "-p",
            "--prefix",
            dest="prefix",
            help="prefix file",
            metavar="FILE")
        (opts, args) = parser.parse_args()
        if opts.config is None:
                parser.print_help()
        else:
                out_prex = opts.prefix

                cf = ConfigParser.ConfigParser()
                cf.read(opts.config)

                SV_config_file = cf.get("pysubsimtree_settings", "SV_config_file")
                ref_fasta = cf.get("pysubsimtree_settings", "ref_fasta")
                dbsnp = cf.get("pysubsimtree_settings", "dbsnp")
                germline_num = cf.getint("pysubsimtree_settings", "germline_num")
                hyp_rate = cf.getfloat("pysubsimtree_settings", 'hyp_rate')
                chrome_can = cf.get("pysubsimtree_settings", "chrome")
                pf_file_name = cf.get("pysubsimtree_settings",
                                      "population_fraction_config_file")
                pf_dict = read_pf_config(pf_file_name)
                ref = reference(ref_fasta)
                snp_dic = read_dbsnp(dbsnp, chrome_can)

                ref_result = generate_normal(
                    copy.deepcopy(ref),
                    snp_dic,
                    germline_num,
                    'germline_SNP_pois.txt',
                    hyp_rate)
                ref_germline = ref_result[0]
                # 生成germline ref
                ref_germline_new = ref_germline
                # 输出germline fasta文件
                # outputFa(ref_germline_new, out_prex+'germline.fa')

                vt = VariantTree(SV_config_file)
                vt.makeSV(ref_germline_new)

                for popu in pf_dict.keys():
                        vt.variant_tree[popu].outputRef(
                            ref_germline_new, out_prex)
                        vt.variant_tree[popu].outputVartPosi(out_prex)


if __name__ == "__main__":
        start = time.clock()
        print('simulation begins at:'+str(start))
        main()
        end = time.clock()
        print('simulation ends at:'+str(end))
        print("The function run time is : %.03f seconds" % (end-start))
