# Pysubsim-tree


## Introduction

Next-generation sequencing (NGS) and the third generation sequencing (TGS) have
recently allowed us to develop algorithms to quantitatively dissect the extent
of het- erogeneity within a tumour, resolve cancer evolution history and
identify the somatic variations and aneuploidy events.  Simulation of tumor NGS
and TGS data serves as a powerful and cost-effective approach for benchmarking
these algorithms, however, there is no available tool that could simulate all
the distinct subclonal genomes with diverse aneuploidy events and somatic
variations according to the given tumor evolution history. We provide a
simulation package, Pysubsim-tree, which could simulate the tumor genomes
according to their evolution history defined by the somatic variations and
aneuploidy events. 


Pysubsim-tree could be run individually or strung with art_illumina, bwa and
samtools to generate Next-Generation Sequencing (NGS) reads according to the
different cancer genomes, which can be mixed with specified proportions.

## Example

	python $PYSUBSIMTREEPATH/pysubsimtree.py -c $CONFIG_PATH/ini.config



