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
from genome_range import GenomeRange

from utils import outputFa

from anytree import NodeMixin
import copy


class VariantNode(NodeMixin):
    def __init__(self, name, parent=None):
        self.name = name
        self.parent = parent
        self.sv_list = []

        # 变异位置
        self.sv_positions = None
        # 断点位置和信息
        self.breakpoints = None
        # 可用位置
        self.avail_position = None

    def outputRef(self, ref, outfilePrefix):
        variant_ref = self.breakpoints.generateFa(ref)
        outputFa(variant_ref, outfilePrefix+self.name+".fa")

    def outputVartPosi(self, outfilePrefix):
        # self.sv_positions[chroms[k]].append([
        # pois,
        # self.sv_list[i][0],
        # self.sv_list[i][1],
        # self.sv_list[i][2],
        # self.sv_list[i][-1]])
        # temp_Node.sv_list.append(
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
        """生成变异
        :returns: TODO

        """
        # 获取父节点中可用空间
        self.avail_position = self._getPrtAvailPst(ref)
        # 获取父节点中变异位置
        self.sv_positions = self._getPrtSVPosis()
        # 获取父节点中breakpoint位置
        self.breakpoints = self._getPrtBreakpoints()

        # 生成当前节点的变异位置，将所有变异位置写入self.sv_positions，返回当前节点的变异位置
        new_positions = self._generateNewPosition(ref)

        self.breakpoints.generateBPs(new_positions, self.avail_position,
                                     ref)

    def _generateNewPosition(self, ref):
        current_sv_positions = {}
        chroms = ref.keys()

        # str_len=[len(ref[tmp][0]) for tmp in dic]
        for i in range(len(self.sv_list)):
            # 此处sv变异被按顺序分配到染色体上
            k = i % len(chroms)
            while True:
                count = 1
                pois = self.avail_position.sample1pois(chroms[k])
                # 生成位置，以随机生成的第一个首字母为起始位置。
                if self.avail_position.isOverlaped(chroms[k], pois,
                                                   pois+self.sv_list[i][0]):
                    if chroms[k] not in self.sv_positions.keys():
                        self.sv_positions[chroms[k]] = []
                    self.sv_positions[chroms[k]].append([
                        pois,
                        self.sv_list[i][0],
                        self.sv_list[i][1],
                        self.sv_list[i][2],
                        self.sv_list[i][-1]])
                    if chroms[k] not in current_sv_positions.keys():
                        current_sv_positions[chroms[k]] = []
                    current_sv_positions[chroms[k]].append([
                        pois,
                        self.sv_list[i][0],
                        self.sv_list[i][1],
                        self.sv_list[i][2],
                        self.sv_list[i][-1]])
                    self.avail_position.takePois(chroms[k], pois,
                                                 pois+self.sv_list[i][0])
                    break
                else:
                    count = count+1
                    if count < 20:
                        continue
                    else:
                        break

        for chrom in self.sv_positions.keys():
            self.sv_positions[chrom] = sorted(self.sv_positions[chrom],
                                              key=lambda d: d[0], reverse=True)
        for chrom in current_sv_positions.keys():
            current_sv_positions[chrom] = sorted(current_sv_positions[chrom],
                                                 key=lambda d: d[0],
                                                 reverse=True)

        return current_sv_positions

    def _getPrtBreakpoints(self):
        if self._isRoot():
            return BreakPoints()
        else:
            return copy.deepcopy(self.parent.breakpoints)
        pass

    def _getPrtSVPosis(self):
        if self._isRoot():
            return {}
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
                variant_name = list_line[1]
                variant_length = int(list_line[2])
                variant_copy_number = int(list_line[3])
                variant_genotype = list_line[4]
                variant_number = int(list_line[5])

                if subclonal_name not in variant_nodes.keys():
                    # 此处添加ref
                    temp_Node = VariantNode(subclonal_name)
                    variant_nodes[subclonal_name] = temp_Node
                else:
                    temp_Node = variant_nodes[subclonal_name]

                for i in range(variant_number):
                    # value: 变异长度、拷贝数、基因型、名称
                    temp_Node.sv_list.append(
                        [variant_length, variant_copy_number,
                         variant_genotype, variant_name])

        # 此处添加根节点，也就是germline
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

        return variant_nodes

    def makeSV(self, ref):
        """生成变异位置
        :returns: TODO

        """
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
    #     print node.name, node.sv_list


if __name__ == "__main__":
    main()
