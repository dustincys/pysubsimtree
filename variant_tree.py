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
from genome_range import GenomeRange

from utils.utils import outputFa

from anytree import NodeMixin

import copy
from collection import Counter


variantNodeType = ["SV", "SNV", "PLOIDY"]


class VariantNode(NodeMixin):
    def __init__(self, name, ref, ploidy_type, parent=None, variant_type="SV"):
        assert variant_type in variantNodeType

        self.name = name
        self.parent = parent
        # 用来传递生成ＳＮＶ
        self.ref = ref

        self.variant_type = variant_type

        self.variant_list = []

        # generate ploidy reference
        self.ploidy_status = {}

        self.sv_positions = None
        self.breakpoints = None
        self.avail_position = None

    def outputRef(self, ref, outfilePrefix):
        variant_ref = self.breakpoints.generateFa(ref)
        outputFa(variant_ref, outfilePrefix+self.name+".fa")

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
        outfileHead =\
            "chrom\tpois\tlength\tcopy_number\tgenotype\tvariant_name\n"

        with open(outfileName, 'w') as outfile:
            outfile.write(outfileHead)
            for chrom in self.sv_positions.keys():
                sv_items = self.sv_positions[chrom]
                for sv_item in sv_items:
                    outfile.write(
                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrom,
                                                                sv_item[0],
                                                                sv_item[1],
                                                                sv_item[2],
                                                                sv_item[3],
                                                                sv_item[4]))
        pass

    def make(self, ref):
        """generate variats
        :returns: TODO

        """
        # 首先初始化节点, 即从父节点继承信息
        self._init_node()

        # 第一步，生成ploidy 变异
        if self.variant_type == "PLOIDY":
            self._make_ploidy()
        else:

            # 第一步，生成非结构变异，是否与结构变异出现重合现象？
            # 生成的位置包括： available中位置， insertion中的位置
            self._make_snv()

        # 第二步，生成结构变异
            self._make_sv()

        new_positions = self._generateNewPosition(ref)
        self.breakpoints.generateBPs(new_positions, self.avail_position,
                                     ref)

    def _init_node(self, ref):
        self.avail_position = self._getPrtAvailPst(ref)
        self.sv_positions = self._getPrtSVPosis()
        self.breakpoints = self._getPrtBreakpoints()

    def _make_ploidy(self):
        # 此处对从父节点继承而来的变异信息进行染色体复制操作
        # 单体型需要用字典来实现ploidy变异
        # 从父节点的breakpoints和sv_positions中生成
        for variant in self.variant_list:
            chrom = variant[0]
            # ploidy_number_before = variant[1]
            # ploidy_number_after = variant[2]
            ploidy_type_before = variant[3]
            ploidy_type_after = variant[4]
            ptc_before = Counter(ploidy_type_before)
            ptc_after = Counter(ploidy_type_before)

            self.ploidy_status[chrom] = ploidy_type_after

            for hapl in ptc_before.keys():
                if hapl in ptc_after.keys():
                    if ptc_before[hapl] > ptc_after[hapl]:
                        # 此处可能不需要针对sv_position进行ploidy
                        # 操作，而只需要对breakpoint 进行操作
                        self.breakpoints.delete_ploidy(
                            chrom, hapl, ptc_before[hapl]-ptc_after[hapl])
                    elif ptc_before[hapl] < ptc_after[hapl]:
                        self.breakpoints.delete_ploidy(
                            chrom, hapl, ptc_after[hapl]-ptc_before[hapl])
                    else:
                        continue
                else:
                    self.breakpoints.delete_ploidy(chrom, hapl, ptc_after[hapl])
        pass

    def _make_snp(self):
        pass

    def _make_sv(self):
        pass

    def _generateNewPosition(self):
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
                    end = posi + 1
                else:
                    end = posi + length

                isOverlap = self.avail_position.isOverlaped(chrom, posi, end)
                if isOverlap:
                    # 此处需要在conf 中给定chrom, 因为不同的染色体的ploidy
                    # 不同，对应的变异的基因型也不同
                    self.avail_position.takePosi(chrom, posi, end)

                    if variant_name == "CNV":
                        copy_number = sv[2]
                        genotype = sv[3]
                        self.sv_positions.add_posi_CNV(
                            chrom, posi, length, copy_number, genotype)
                        current_sv_positions.add_posi_CNV(
                            chrom, posi, length, copy_number, genotype)

                    if variant_name == "INVERSION":
                        genotype = sv[2]
                        self.sv_positions.add_posi_INVERSION(
                            chrom, posi, length, genotype)
                        current_sv_positions.add_posi_INVERSION(
                            chrom, posi, length, genotype)

                    if variant_name == "DELETION":
                        genotype = sv[2]
                        self.sv_positions.add_posi_DELETION(
                            chrom, posi, length, genotype)
                        current_sv_positions.add_posi_DELETION(
                            chrom, posi, length, genotype)

                    if variant_name == "TRANSVERSION":
                        genotype = sv[2]
                        self.sv_positions.add_posi_TRANSVERSION(
                            chrom, posi, length, genotype)
                        current_sv_positions.add_posi_TRANSVERSION(
                            chrom, posi, length, genotype)
                    break
                else:
                    count = count+1
                    if count < 20:
                        continue
                    else:
                        break

        self.SV_positions.sorted()
        current_sv_positions.sorted()

        return current_sv_positions

    def _getPrtBreakpoints(self):
        if self._isRoot():
            return BreakPoints()
        else:
            return copy.deepcopy(self.parent.breakpoints)
        pass

    def _getPrtSVPosis(self):
        if self._isRoot():
            return SV_positions()
        else:
            return copy.deepcopy(self.parent.sv_positions)
        pass

    def _getPrtAvailPst(self, ref):
        if self._isRoot():
            return GenomeRange(ref)
        else:
            return copy.deepcopy(self.parent.avail_position)

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
                        variant_length = int(list_line[4])
                        variant_genotype = list_line[5]
                        number = int(list_line[6])

                        variant = [
                            chrom,
                            variant_length,
                            variant_genotype,
                            variant_name,
                            variant_type]
                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)

                    elif variant_name == "INSERTION":
                        chrom = list_line[3]
                        variant_length = int(list_line[4])
                        variant_genotype = list_line[5]
                        ploidy_genotype = list_line[6]
                        number = int(list_line[7])

                        variant = [
                            chrom,
                            variant_length,
                            variant_genotype,
                            ploidy_genotype,
                            variant_name,
                            variant_type]
                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)

                    elif variant_name == "DELETION":
                        chrom = list_line[3]
                        variant_length = int(list_line[4])
                        variant_genotype = list_line[5]
                        ploidy_genotype = list_line[6]
                        number = int(list_line[7])

                        variant = [
                            chrom,
                            variant_length,
                            variant_genotype,
                            ploidy_genotype,
                            variant_name,
                            variant_type]
                        self._add2node(
                            variant, number, subclonal_name, variant_nodes)
                    elif variant_name == "TRANSVERSION":
                        chrom = list_line[3]
                        variant_length = int(list_line[4])
                        variant_genotype = list_line[5]
                        ploidy_genotype = list_line[6]
                        number = int(list_line[7])

                        variant = [
                            variant_length,
                            variant_genotype,
                            ploidy_genotype,
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
                               ploidy_type_after]

                    self._add2node(
                        variant, 1, subclonal_name, variant_nodes)

        self._linkNode(variant_nodes)

        return variant_nodes

    def _add2node(self, variant, number, subclonal_name, variant_nodes):
        if subclonal_name not in variant_nodes.keys():
            temp_Node = VariantNode(subclonal_name)
            variant_nodes[subclonal_name] = temp_Node
        else:
            temp_Node = variant_nodes[subclonal_name]

        for i in range(number):
            temp_Node.variant_list = temp_Node.variant_list + variant

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
