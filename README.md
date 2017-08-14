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

## Usage

	python $PYSUBSIMTREEPATH/pysubsimtree.py -c $CONFIG_PATH/ini.config -p $PREFIX

## Configuration

A full demo is given in `./test`.
In Pysubsim-tree, Node labelled `1` is normal genome.
Node labelled `1.1` is tumor subclone genome, Node `1` is  parent node of Node `1.1` Node `1.1`. 


`./test/config/init.conf` 	main configuration file


	[pysubsimtree_settings]
	#	Tumor evolution tree configuration file
	SV_config_file=/media/d/github/pysubsimtree/test/config/variantTree.conf
	#	reference genome
	ref_fasta=/media/d/github/pysubsimtree/test/small.fa
	#	dbsnp file
	dbsnp=/media/d/github/pysubsimtree/test/small.vcf
	#	number of germline snps
	germline_num=2
	#	ratio of heterozygou germline snps
	hyp_rate=0.5
	#	chromosome
	chrome=chr1
	# 	Tumor subpopulation fraction configuration file
	population_fraction_config_file=/media/d/github/pysubsimtree/test/config/pf_config


`./test/config/variantTree.conf` 	Tumor evolution tree configuration file


	# This is variant tree configuration file
	# The sv configuration
	#
	# CNV
	# Variant node label	Variant node type	Variant type	chrom	Length of CNVs	Copy number	Genotype	Number
	1.1	SV	CNV	chr1	2	0	NONE	1
	1.3.1	SV	CNV	chr1	2	0	NONE	1
	1.3.1	SV	CNV	chr1	2	1	P	1
	#1.1	SV	CNV	chr1	2	1	M	1
	1.3.1	SV	CNV	chr1	2	2	PP	1
	1.2	SV	CNV	chr1	2	2	MM	1
	1.2.1	SV	CNV	chr1	2	3	PPM	1
	#1.2.1	SV	CNV	chr1	2	3	PMM	1
	#1.2.1	SV	CNV	chr1	2	3	MMM	1
	#1.2.1	SV	CNV	chr1	2	4	PPPP	1
	#1.2.1	SV	CNV	chr1	2	4	PPPM	1
	#1.3	SV	CNV	chr1	2	4	PPMM	1
	#1.3	SV	CNV	chr1	2	4	PMMM	1
	1.3.1	SV	CNV	chr1	2	4	MMMM	1
	#
	# INVERTION
	# Variant node label	Variant node type	Variant type	chrom	haplotype(Paternal or Maternal)	index of chromatid	Length	Genotype	Number
	1.3.1	SV	INVERTION	chr1	P	0	4	2
	#
	# TANDEMDUP
	# Variant node label	Variant node type	Variant type	chrom	haplotype(Paternal or Maternal)	index of chromatid	Length	Number
	1.3.1	SV	TANDEMDUP	chr1	P	0	4	2	2
	#
	# INSERTION
	# Variant node label	Variant node type	Variant type	chrom	haplotype(Paternal or Maternal)	index of chromatid	Length	Genotype	Number
	1.3.1	SV	INSERTION	chr1	P	0	4	2
	#
	# DELETION
	# Variant node label	Variant node type	Variant type	chrom	haplotype(Paternal or Maternal)	index of chromatid	Length	Genotype	Number
	1.3.1	SV	DELETION	chr1	M	0	4	2
	#
	# TRANSLOCATION
	# Variant node label	Variant node type	Variant type	chrom from	haplotype(Paternal or Maternal) from	index of chromatid from	Length	chrom to	haplotype(Paternal or Maternal) to	index of chromatid to	Number
	1.3.1	SV	TRANSLOCATION	chr1	M	0	4	chr1	P	0	2
	#
	# SNV
	# Variant node label	Variant node type	chrom	is heterozygous	is overlaped	number
	1.3.1	SNV	chr1	FALSE	TRUE	2	2
	1.3.1	SNV	chr1	TRUE	TRUE	2	2
	1.3.1	SNV	chr1	FALSE	FALSE	2	2
	#
	# Aneuploidy event
	# Variant node label	Variant node type	chrom 	ploidy number before	ploidy number after	ploidy status before	ploidy status after
	1.3	PLOIDY	chr1	2	1	PM	PPM


`./test/config/pf_config` 	Tumor subpopulation fraction configuration file

	#name	fraction
	1	0.2
	1.1	0.3
	1.2	0.2
	1.2.1	0.1
	1.3	0.1
	1.3.1	0.1
