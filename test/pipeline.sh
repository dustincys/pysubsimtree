#!/bin/sh

PYSUBSIMTREEPATH=/media/d/github/pysubsimtree
REFERENCE_PATH=/media/d/github/pysubsimtree

CONFIG_PATH=./config


################  specify_ploidy  #########################################
#python $PYSUBSIMTREEPATH/specify_ploidy.py -i $REFERENCE_PATH/chr1.fa  -c $CONFIG_PATH/ploidy.conf -o output_ploidy.fa
#
##################  generate sv fasta  #########################################

python $PYSUBSIMTREEPATH/pysubsimtree.py -c $CONFIG_PATH/init.config

#
##################  generate fastq  #########################################
#
#python  $PYSUBSIMTREEPATH/run_art.py -i ./test_chr1_germline.fa -o germline
#python  $PYSUBSIMTREEPATH/run_art.py -i ./test_chr1_somatic.fa -o subclone1
#python  $PYSUBSIMTREEPATH/run_art.py -i ./test_chr1_somatic_subclone_1.fa -o subclone2
#python  $PYSUBSIMTREEPATH/run_art.py -i ./test_chr1_somatic_subclone_2.fa -o subclone3

################  generate mixture fq  #########################################

#python $PYSUBSIMTREEPATH/mixture_v0.5.py -i ./config/subclonepurity.conf -o mixture_tumor

# bwa mem $REFERENCE_PATH/chr1.fa mixture_tumor_1.fq mixture_tumor_2.fq -t 30  > mixture_tumor.sam
# samtools view -Sb mixture_tumor.sam > mixture_tumor.bam
# rm mixture_tumor.sam
# samtools sort mixture_tumor.bam -o mixture_tumor.sorted.bam
# samtools index mixture_tumor.sorted.bam
