import sys
from optparse import OptionParser
import getopt
import time
import os
def read_config(config_filename):
    dic={}
    for line in open(config_filename):
        if not line.startswith('#'):
            newline=line.rstrip().split('\t')
            dic[newline[0]]=int(newline[1])
        else:
            continue
    return dic
def ploidy(refDic,dic):
    refDicNew={}
    for key in dic:
        for n in range(dic[key]):
            refDicNew[key+'-hap-'+str(n+1)]=refDic[key]
    return refDicNew
def read_fasta(filename):
    refDic={}
    chrName=''
    for line in open(filename):
        newline=line.rstrip()
        if newline.startswith('>'):
            if chrName!='':
                if not chrName.startswith('chr'):
                    chrName='chr'+chrName
                refDic[chrName]=tmpStr
            chrName=newline.split('>')[1]
            if not chrName.startswith('chr'):
                chrName='chr'+chrName
            tmpStr=''
        else:
            tmpStr=tmpStr+newline.upper()
    refDic[chrName]=tmpStr
    return refDic
def output(ref,outfilename):
    outfile=open(outfilename,'w')
    tmpKey=sorted(ref.keys())
    for key in tmpKey:
        i=0
        strLen=len(ref[key])
        outfile.write('>'+key+'\n')
        while i+50<=strLen:
            outfile.write(ref[key][i:i+50]+'\n')
            i=i+50
    outfile.close()
def main():
    usage = """%prog -i <file> -c <config file>  -o <out_fasta>

specify_ploidy
Author: Yuchao Xia
Description: specify the number of ploidy for different chromesome
	"""

    parser = OptionParser(usage)
    parser.add_option("-i", "--inFile", dest="inFile", help="A reference fasta file.",metavar="FILE")
    parser.add_option("-c","--config",dest='config',help='the number of ploidy of different chromesome',metavar='FILE')
    parser.add_option('-o','--output',dest='output',help='output fasta file',metavar="file")

    (opts, args) = parser.parse_args()
    if opts.inFile is None or opts.config is None:
        parser.print_help()
    else:
        if opts.output is None:
            outfilename='output_ploid.fa'
        else:
            outfilename=opts.output
        config_dic=read_config(opts.config)
        refDic=read_fasta(opts.inFile)
        refDic=ploidy(refDic,config_dic)
        output(refDic,outfilename)
if __name__ == "__main__":
    print time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    start = time.clock()
    print('simulation begins at:'+str(start))
    main()
    end = time.clock()
    print('simulation ends at:'+str(end))
    print("The function run time is : %.03f seconds" %(end-start))
    print time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
