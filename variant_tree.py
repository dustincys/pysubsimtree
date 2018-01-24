#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: variant_tree.py
#          Desc: 根据配置文件，生成变异树，根据要求可以遍历树的每一个节点
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-07-25 22:45:00
#       History:
# =============================================================================
'''

from breakpoint import BreakPoints
from sv_positions import SV_positions
from snv_positions import SNV_positions
from genome_range import GenomeRange

from utils.utils import outputFa, SNP
from utils import constants

from anytree import NodeMixin

import copy
from collections import Counter
import random


variantNodeType = ["SV", "SNV", "PLOIDY"]


class VariantNode(NodeMixin):
    def __init__(self, name, parent=None, variant_type="SV"):
        assert variant_type in variantNodeType

        self.name = name
        self.parent = parent
        # 用来传递生成ＳＮＶ

        self.variant_type = variant_type

        self.variant_list = []

        # generate ploidy reference
        self.ploidy_status = {}

        self.snv_positions = None

        self.sv_positions = None
        self.breakpoints = None
        self.avail_position = None

        # record the none del sv position
        self.noDEL_position = None

    def outputRef(self, ref, outfilePrefix):
        ref_ploidy = self._ref_make_ploidy(copy.deepcopy(ref))
        ref_ploidy_snv = self._ref_make_snv(ref_ploidy)
        variant_ref = self.breakpoints.generateFa(ref_ploidy_snv)
        outputFa(variant_ref, outfilePrefix+self.name+".fa")

    def _ref_make_ploidy(self, ref):
        ref_ploidy = {}
        for chrom in self.ploidy_status.keys():
            if chrom not in ref_ploidy.keys():
                ref_ploidy[chrom] = {}

            psc = Counter(self.ploidy_status[chrom])
            for hapl_type in psc.keys():
                if hapl_type not in ref_ploidy[chrom].keys():
                    ref_ploidy[chrom][hapl_type] = []

                for i in range(psc[hapl_type]):
                    ref_ploidy[chrom][hapl_type].append(
                        ref[chrom][hapl_type][0])

        return ref_ploidy

    def _ref_make_snv(self, ref_ploidy):
        # 与breakpoint结构一致
        for chrom in self.snv_positions.snvp_dict.keys():
            for hapl_type in self.snv_positions.snvp_dict[chrom].keys():
                for hapl_idx in \
                        range(len(self.snv_positions.snvp_dict
                                  [chrom][hapl_type])):
                    lstr = list(ref_ploidy[chrom][hapl_type][hapl_idx])
                    for snvp in\
                            self.snv_positions.snvp_dict[
                                                chrom][hapl_type][hapl_idx]:
                        lstr[snvp.position] = snvp.B_allele

                    ref_ploidy[chrom][hapl_type][hapl_idx] = ''.join(lstr)

        return ref_ploidy

    def outputVartPosi(self, outfilePrefix):
        # self.sv_positions[chroms[k]].append([
        # posi,
        # self.variant_list[i][0],
        # self.variant_list[i][1],
        # self.variant_list[i][2],
        # self.variant_list[i][-1]])
        # temp_Node.variant_list.append(
        # [variant_length, variant_copy_number,
        # variant_genotype, variant_name])
        outfileName = outfilePrefix+self.name+".SV.txt"

        with open(outfileName, 'w') as outfile:
            for chrom in self.sv_positions.svp_dict.keys():
                sv_items = self.sv_positions.svp_dict[chrom]
                for sv_item in sv_items:
                    outfile.write("chrom\t" + sv_item.info_str_title())
                    outfile.write(chrom + "\t" + sv_item.info_str())
        pass

    def make(self, ref):
        """generate variats
        :returns: TODO

        """
        # 首先初始化节点, 即从父节点继承信息
        self._init_node(ref)

        # 第一步，生成ploidy 变异
        if self.variant_type == "PLOIDY":
            self._make_ploidy()
        else:

            # 第一步，生成非结构变异，是否与结构变异出现重合现象？
            # 生成的位置包括： available中位置， insertion中的位置

            ref_ploidy = self._ref_make_ploidy(copy.deepcopy(ref))
            ref_ploidy_snv = self._ref_make_snv(ref_ploidy)

            # 第二步，生成结构变异, 此处是生成了ploidy的ref
            self._make_sv(ref_ploidy_snv)
            self._make_snv(ref_ploidy_snv)

    def _init_node(self, ref):
        self.avail_position = self._getPrtAvailPst(ref)
        self.sv_positions = self._getPrtSVPosis()
        self.snv_positions = self._getPrtSNVPosis()
        self.breakpoints = self._getPrtBreakpoints()
        self.ploidy_status = self._getPrtPloidyStatus(ref)
        self.noDEL_position = self._getPrtNBPPst(ref)

    def _make_ploidy(self):
        # 此处对从父节点继承而来的变异信息进行染色体复制操作
        # 单体型需要用字典来实现ploidy变异
        # 从父节点的breakpoints和sv_positions中生成

        # variant = [chrom, ploidy_number_before,
                    # ploidy_number_after, ploidy_type_before,
                    # ploidy_type_after]

        for variant in self.variant_list:
            chrom = variant[0]
            # ploidy_number_before = variant[1]
            # ploidy_number_after = variant[2]
            ploidy_type_before = variant[3]
            ploidy_type_after = variant[4]
            ptc_before = Counter(ploidy_type_before)
            ptc_after = Counter(ploidy_type_after)

            for hapl in ptc_before.keys():
                if hapl in ptc_after.keys():
                    if ptc_before[hapl] > ptc_after[hapl]:
                        # 此处可能不需要针对sv_position进行ploidy
                        # 操作，而只需要对breakpoint 进行操作
                        number = ptc_before[hapl]-ptc_after[hapl]
                        hapl_idxes = random.sample(
                            range(ptc_before[hapl]), number)
                        self.breakpoints.delete_ploidy(chrom, hapl, hapl_idxes)
                        self.snv_positions.delete_ploidy(chrom, hapl,
                                                         hapl_idxes)

                    elif ptc_before[hapl] < ptc_after[hapl]:
                        number = ptc_after[hapl]-ptc_before[hapl]
                        hapls = [random.choice(range(ptc_before[hapl])) for _ in
                                 range(number)]
                        self.breakpoints.add_ploidy(chrom, hapl, hapls)
                        self.snv_positions.add_ploidy(chrom, hapl, hapls)

                    else:
                        continue
                else:
                    number = ptc_before[hapl]
                    hapl_idxes = range(ptc_before[hapl])
                    self.breakpoints.delete_ploidy(chrom, hapl, hapl_idxes)
                    self.snv_positions.delete_ploidy(chrom, hapl, hapl_idxes)

        self.ploidy_status[chrom] = ploidy_type_after
        pass

    def _make_snv(self, ref):
        # 从avail_position 中生成非重合变异，
        # 从sv_position和breakpoints的insertStr中生成重合的变异
        # 只有非重合变异可以生成纯合突变

        snv_list = filter(lambda item: item[-1] == "SNV", self.variant_list)
        for snv in snv_list:
            chrom = snv[0]
            isHetero = snv[1]
            isOverlap = snv[2]
            number = snv[3]

            if isOverlap == "TRUE":
                # only hetero
                number_noDEL = number
                number_bp = 0
                noDELCNVpois = self._getNDCPs(chrom)
                strBps = self._getStrBps(chrom)

                if strBps is not None and len(noDELCNVpois) > 0:
                    number_noDEL = int(number / 2)
                    number_bp = number - number_noDEL
                if len(noDELCNVpois) == 0:
                    number_noDEL = 0
                    number_bp = number
                if strBps is None:
                    number_bp = 0
                    number_noDEL = number
                if strBps is None and len(noDELCNVpois) == 0:
                    number_bp = 0
                    number_noDEL = 0

                for i in range(number_bp):
                    bp = random.sample(strBps, 1)[0]
                    position = random.sample(range(len(bp.insertStr)), 1)[0]
                    bp_lstr = list(bp.insertStr)
                    bp_lstr[position] = SNP(bp_lstr[position])
                    bp.insertStr = ''.join(bp_lstr)

                for i in range(number_noDEL):

                    noDELCNV = random.sample(noDELCNVpois, 1)[0]
                    length = noDELCNV.sv.length
                    hapl_remain = noDELCNV.sv.hapl_remain
                    position = noDELCNV.position\
                        + random.sample(range(length), 1)[0]

                    temp_hapl_type = random.sample(hapl_remain.keys(), 1)[0]
                    temp_hapl_index = random.sample(
                        hapl_remain[temp_hapl_type], 1)[0]
                    hapl_lstr = ref[chrom][temp_hapl_type][temp_hapl_index]
                    B_allele = SNP(hapl_lstr[position])

                    self.snv_positions.addPosi(chrom, temp_hapl_type,
                                               temp_hapl_index,
                                               position, B_allele, isHetero,
                                               isOverlap)
            else:

                positions = self.avail_position.samplePosis(chrom, number)

                if isHetero == "TRUE":
                    hapl_type = random.sample(ref[chrom].keys(), 1)[0]

                    for pi in positions:
                        hapl_index = random.sample(
                            range(len(ref[chrom][hapl_type])), 1)[0]
                        hapl_lstr = ref[chrom][hapl_type][hapl_index]
                        B_allele = SNP(hapl_lstr[pi])

                        self.snv_positions.addPosi(chrom, hapl_type, hapl_index,
                                                   pi, B_allele, isHetero,
                                                   isOverlap)
                else:
                    for pi in positions:
                        temp_hapl_type = random.sample(ref[chrom].keys(), 1)[0]
                        temp_hapl_index = random.sample(
                            range(len(ref[chrom][temp_hapl_type])), 1)[0]
                        B_allele = SNP(
                            ref[chrom][temp_hapl_type][temp_hapl_index][pi]
                        )
                        for hapl_type in ref[chrom].keys():
                            for hapl_index in range(len(ref[chrom][hapl_type])):
                                self.snv_positions.addPosi(chrom, hapl_type,
                                                           hapl_index,
                                                           pi, B_allele,
                                                           isHetero, isOverlap)
        pass

    def _getStrBps(self, chrom):
        # 此处为了进行ploidy操作，bp结构为 {chr1:{'P':[[bp1],[bp2]], 'M':[[bp3]]}}

        for hapl_type in self.breakpoints.breaks_list[chrom].keys():
            for hapl_index in range(len(
                    self.breakpoints.breaks_list[chrom][hapl_type])):
                strBps = filter(lambda item: item.insertStr is not None,
                                self.breakpoints.breaks_list
                                [chrom][hapl_type][hapl_index])
                if len(strBps) > 1:
                    return strBps

        return None

    def _getNDCPs(self, chrom):
        return filter(
            lambda item: item.sv_type == "CNV" and item.sv.genotype != "NONE",
            self.sv_positions.svp_dict[chrom])

    def _make_sv(self, ref):
        new_positions = self._generateNewPosition()
        # noBP remove deletion seg
        self.breakpoints.generateBPs(new_positions, self.avail_position,
                                     self.noDEL_position, ref,
                                     self.ploidy_status)

    def _generateNewPosition(self):
        clmin = constants.COMPLEXINDEL_LENGTH_MIN
        clmax = constants.COMPLEXINDEL_LENGTH_MAX

        current_sv_positions = SV_positions()

        # str_len=[len(ref[tmp][0]) for tmp in dic]
        sv_list = filter(lambda item: item[-1] == "SV", self.variant_list)

        for sv in sv_list:
            # 此处sv变异被按顺序分配到染色体上
            chrom = sv[0]
            length = sv[1]
            variant_name = sv[-2]

            while True:
                count = 1
                posi = self.avail_position.sample1posi(chrom)
                # 生成位置，以随机生成的第一个首字母为起始位置。
                if variant_name == "INSERTION":
                    isOverlap = self.avail_position.isOverlaped(chrom, posi,
                                                                posi + 1)
                else:
                    isOverlap = self.avail_position.isOverlaped(chrom, posi,
                                                                posi + length)

                if isOverlap:
                    # 此处需要在conf 中给定chrom, 因为不同的染色体的ploidy
                    # 不同，对应的变异的基因型也不同
                    if variant_name == "INSERTION":
                        self.avail_position.takePosi(chrom, posi, posi + 1)
                    else:
                        self.avail_position.takePosi(chrom, posi, posi + length)

                    if variant_name == "CNV":
                        copy_number = sv[2]
                        genotype = sv[3]

                        self.sv_positions.add_posi_CNV(
                            chrom, posi, length, copy_number, genotype,
                            self.ploidy_status)
                        current_sv_positions.add_posi_CNV(
                            chrom, posi, length, copy_number, genotype,
                            self.ploidy_status)

                    if variant_name == "INVERTION":

                        hapl_type = sv[2]
                        hapl_idx = sv[3]
                        # ploidy_status = sv[4]

                        self.sv_positions.add_posi_INVERTION(
                            chrom, hapl_type, hapl_idx, posi, length)
                        current_sv_positions.add_posi_INVERTION(
                            chrom, hapl_type, hapl_idx, posi, length)

                    if variant_name == "TANDEMDUP":

                        hapl_type = sv[2]
                        hapl_idx = sv[3]
                        times = sv[4]
                        # ploidy_status = sv[5]

                        self.sv_positions.add_posi_TANDEMDUP(
                            chrom, hapl_type, hapl_idx, posi, length, times)
                        current_sv_positions.add_posi_TANDEMDUP(
                            chrom, hapl_type, hapl_idx, posi, length, times)

                    if variant_name == "DELETION":

                        hapl_type = sv[2]
                        hapl_idx = sv[3]
                        # ploidy_status = sv[4]

                        self.sv_positions.add_posi_DELETION(
                            chrom, hapl_type, hapl_idx, posi, length)
                        current_sv_positions.add_posi_DELETION(
                            chrom, hapl_type, hapl_idx, posi, length)

                    if variant_name == "COMPLEXINDEL":

                        hapl_type = sv[2]
                        hapl_idx = sv[3]
                        # ploidy_status = sv[4]
                        length1 = random.randint(clmin, clmax)
                        length2 = random.randint(clmin, clmax)

                        self.sv_positions.add_posi_COMPLEXINDEL(
                            chrom, hapl_type, hapl_idx, posi, length1, length2)
                        current_sv_positions.add_posi_COMPLEXINDEL(
                            chrom, hapl_type, hapl_idx, posi, length1, length2)

                    if variant_name == "INSERTION":

                        hapl_type = sv[2]
                        hapl_idx = sv[3]
                        # ploidy_status = sv[4]

                        self.sv_positions.add_posi_INSERTION(
                            chrom, hapl_type, hapl_idx, posi, length)
                        current_sv_positions.add_posi_INSERTION(
                            chrom, hapl_type, hapl_idx, posi, length)

                    if variant_name == "TRANSLOCATION":
                        hapl_type_from = sv[2]
                        hapl_idx_from = int(sv[3])
                        chrom_to = sv[4]
                        hapl_type_to = sv[5]
                        hapl_idx_to = int(sv[6])
                        # ploidy_genotype = sv[7]

                        self.sv_positions.add_posi_TRANSLOCATION(
                            chrom, posi, hapl_type_from, hapl_idx_from,
                            chrom_to, hapl_type_to, hapl_idx_to, length)

                        current_sv_positions.add_posi_TRANSLOCATION(
                            chrom, posi, hapl_type_from, hapl_idx_from,
                            chrom_to, hapl_type_to, hapl_idx_to, length)
                    break
                else:
                    count = count+1
                    if count < 20:
                        continue
                    else:
                        break

        self.sv_positions.sorted()
        current_sv_positions.sorted()

        return current_sv_positions

    def _getPrtBreakpoints(self):
        if self._isRoot():
            return BreakPoints()
        else:
            return copy.deepcopy(self.parent.breakpoints)
        pass

    def _getPrtSNVPosis(self):
        if self._isRoot():
            return SNV_positions()
        else:
            return copy.deepcopy(self.parent.snv_positions)
        pass

    def _getPrtSVPosis(self):
        if self._isRoot():
            return SV_positions()
        else:
            return copy.deepcopy(self.parent.sv_positions)
        pass

    def _getPrtAvailPst(self, ref):
        if self._isRoot():
            avail_position = GenomeRange()
            avail_position.init_from_ref(ref)
            return avail_position
        else:
            return copy.deepcopy(self.parent.avail_position)

    def _getPrtNBPPst(self, ref):
        if self._isRoot():
            avail_position = GenomeRange()
            avail_position.init_none(ref)
            return avail_position
        else:
            return copy.deepcopy(self.parent.noDEL_position)

    def _getPrtPloidyStatus(self, ref):
        if self._isRoot():
            ps = {}
            for chrom in ref.keys():
                ps[chrom] = "PM"
            return ps
        else:
            return copy.deepcopy(self.parent.ploidy_status)

    def _isRoot(self):
        if self.parent is None:
            return True
        else:
            return False


class VariantTree(object):

    """变异树"""

    def __init__(self, configFileName):
        self.variant_tree = self._read_config(configFileName)

    def _getParentName(self, name):
        """生成父节点的名称

        格式：
        1
        1.1
        1.2
        1.1.1
        ...

        :name: 当前节点名称
        :returns: 父节点名称

        """
        name_list = name.rsplit(".", 1)
        if len(name_list) == 1:
            parent_name = None
        else:
            parent_name = name_list[0]

        return parent_name

    def _read_config(self, fileName):
        variant_nodes = {}

        with open(fileName) as infile:
            for line in sorted(infile):
                # 这里按顺序生成
                if not line.startswith('#'):
                    list_line = line.rstrip().split('\t')
                else:
                    continue

                subclonal_name = list_line[0]
                variant_type = list_line[1]

                if variant_type == "SV":
                    variant_name = list_line[2]

                    if variant_name == "CNV":
                        chrom = list_line[3]
                        variant_length = int(list_line[4])
                        variant_copy_number = int(list_line[5])
                        variant_genotype = list_line[6]
                        number = int(list_line[7])
                        variant = [
                            chrom,
                            variant_length,
                            variant_copy_number,
                            variant_genotype,
                            variant_name,
                            variant_type]
                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)

                    elif variant_name == "INVERTION":
                        chrom = list_line[3]
                        hapl_type = list_line[4]
                        hapl_idx = int(list_line[5])
                        variant_length = int(list_line[6])
                        number = int(list_line[7])

                        variant = [
                            chrom,
                            variant_length,
                            hapl_type,
                            hapl_idx,
                            variant_name,
                            variant_type]

                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)

                    elif variant_name == "TANDEMDUP":
                        chrom = list_line[3]
                        hapl_type = list_line[4]
                        hapl_idx = int(list_line[5])
                        variant_length = int(list_line[6])
                        times = int(list_line[7])
                        number = int(list_line[8])

                        variant = [
                            chrom,
                            variant_length,
                            hapl_type,
                            hapl_idx,
                            times,
                            variant_name,
                            variant_type]

                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)

                    elif variant_name == "INSERTION":
                        chrom = list_line[3]
                        hapl_type = list_line[4]
                        hapl_idx = int(list_line[5])
                        variant_length = int(list_line[6])
                        number = int(list_line[7])

                        variant = [
                            chrom,
                            variant_length,
                            hapl_type,
                            hapl_idx,
                            variant_name,
                            variant_type]

                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)

                    elif variant_name == "DELETION":
                        chrom = list_line[3]
                        hapl_type = list_line[4]
                        hapl_idx = int(list_line[5])
                        variant_length = int(list_line[6])
                        number = int(list_line[7])

                        variant = [
                            chrom,
                            variant_length,
                            hapl_type,
                            hapl_idx,
                            variant_name,
                            variant_type]

                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)

                    elif variant_name == "COMPLEXINDEL":
                        chrom = list_line[3]
                        hapl_type = list_line[4]
                        hapl_idx = int(list_line[5])
                        number = int(list_line[6])

                        variant = [
                            chrom,
                            -1,
                            hapl_type,
                            hapl_idx,
                            variant_name,
                            variant_type]

                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)

                    elif variant_name == "TRANSLOCATION":
                        chrom_from = list_line[3]
                        hapl_type_from = list_line[4]
                        hapl_idx_from = int(list_line[5])
                        variant_length = int(list_line[6])
                        chrom_to = list_line[7]
                        hapl_type_to = list_line[8]
                        hapl_idx_to = int(list_line[9])
                        number = int(list_line[10])

                        variant = [
                            chrom_from,
                            variant_length,
                            hapl_type_from,
                            hapl_idx_from,
                            chrom_to,
                            hapl_type_to,
                            hapl_idx_to,
                            variant_name,
                            variant_type]

                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)
                    else:
                        pass
                elif variant_type == "PLOIDY":
                    chrom = list_line[2]
                    ploidy_number_before = int(list_line[3])
                    ploidy_number_after = int(list_line[4])
                    ploidy_type_before = list_line[5]
                    ploidy_type_after = list_line[6]
                    variant = [chrom, ploidy_number_before,
                               ploidy_number_after, ploidy_type_before,
                               ploidy_type_after, variant_type]

                    self._add2node(
                        variant, 1, subclonal_name, variant_nodes)

                elif variant_type == "SNV":
                    chrom = list_line[2]
                    isHetero = list_line[3]
                    isOverlap = list_line[4]
                    number = int(list_line[6])
                    variant = [chrom, isHetero, isOverlap, number,
                               variant_type]
                    self._add2node(
                        variant, 1, subclonal_name, variant_nodes)

        self._linkNode(variant_nodes)

        return variant_nodes

    def _add2node(self, variant, number, subclonal_name, variant_nodes):
        if subclonal_name not in variant_nodes.keys():
            temp_Node = VariantNode(subclonal_name,
                                    variant_type=variant[-1])
            variant_nodes[subclonal_name] = temp_Node
        else:
            temp_Node = variant_nodes[subclonal_name]

        for i in range(number):
            temp_Node.variant_list = temp_Node.variant_list + [variant]

    def _linkNode(self, variant_nodes):
        if "1" not in variant_nodes.keys():
            variant_nodes["1"] = VariantNode("1")

        for node_name in variant_nodes.keys():
            node_parent_name = self._getParentName(node_name)
            if node_parent_name is not None:
                try:
                    variant_nodes[node_name].parent =\
                        variant_nodes[node_parent_name]
                except KeyError:
                    print "Node parent error!"

    def makeSV(self, ref):
        self._make(self.variant_tree['1'], ref)

    def _make(self, node, ref):
        node.make(ref)
        for child in node.descendants:
            self._make(child, ref)


def main():
    """
    :returns: TODO

    """
    # variant_nodes = read_config("../config/sv_config_test")
    # print RenderTree(variant_nodes['1'])
    # for pre, _, node in RenderTree(variant_nodes['1']):
    #     print node.name, node.variant_list


if __name__ == "__main__":
    main()
