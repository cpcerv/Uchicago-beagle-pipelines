#!/bin/bash
import os, datetime, glob, time, argparse, shutil
parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string)
	return os.path.abspath(string)

parser.add_argument("-d", "--bamDir", type=readable_dir,
		help="directory of bam files")
	
args = parser.parse_args()

outDir = args.bamDir + '_bam_summary'+'.txt'


FilesList= glob.glob(os.path.join(args.bamDir, '*.bam'))
for fs in FilesList:
	f=fs.strip().split('/')[-1]
	os.system('echo -n -e ' + f + '"\t" >> ' + outDir+ '\n')
	os.system('samtools view -F 0x904 -c ' + fs + ' >> ' + outDir + ' \n')


