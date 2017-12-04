#!/usr/bin/bash

import os, datetime, glob, time, argparse

parser = argparse.ArgumentParser(add_help=True)
    
def readable_dir(string):
    if not os.path.isdir(string):
        raise argparse.ArgumentTypeError("%s directory not found" % string)
    return os.path.abspath(string)


parser.add_argument("-d", "--htseqDir", type=readable_dir,
        help="directory of htseq count files, will also make summary output there")
                            
args = parser.parse_args()
readsCount = args.htseqDir+'/htseq_summary.txt'

htseqDir = args.htseqDir
        
DirList= glob.glob(os.path.join(htseqDir,'*[!.txt][!log]*'))

reads_mapped = {}
fout = open(readsCount, 'w')
for fs in DirList:
    name = fs.strip().split('/')[-1].strip().split('.')[0]
    if not ("Undetermined" in fs):
        if not (name in reads_mapped):
            reads_mapped[name] = 0
        fin = open(fs)
        for line in fin:
            if not line.startswith("_"):
                reads_mapped[name] = reads_mapped[name] + int(line.strip().split('\t')[-1])
        fin.close()

fout.write("sample_id\thtseq_assigned\n")
for key in reads_mapped:
    fout.write(key + '\t' + str(reads_mapped[key]) + '\n')
fout.close()
