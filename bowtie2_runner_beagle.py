#!/bin/bash

##### WRITE PBS SCRIPT AND CALL bowite2 #########
#	This python script will write a .pbs for every fastq file      #
#	in a directory.  It calls a master tophat-runner shell script #
#	for aprun.  n and N for aprun will always be 1.                 #
#	Open source. By:Carlos PerezCervantes. Moskowitz lab # 
#####################################################

import os, datetime, glob, time, argparse

parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string)
	return os.path.abspath(string)


parser.add_argument("-d", "--fastqDir", type=readable_dir,
		help="directory of Fastq files")
				  
parser.add_argument("-o", "--out", type=str,
		help="specify output destination")
				  
parser.add_argument("-g", "--genome", type=str,default='mm10',
		help="specify genome: mm9, mm10, hg19, hg38")
				  
parser.add_argument("-p", "--processors", type=str, default='32',
		help="specify number of processors for computation node (a mutiple of 32) and tophat")              

parser.add_argument("-w", "--walltime", type=str, default='8', 
		help="specify running time in hours, i.e 02, 4, 16")

parser.add_argument("-m", "--mapping", type=str, default='1', 
		help="specify number of alignment matches ")
		
args = parser.parse_args()
args.out = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)

jobsDir = args.out+'/log/'

if not os.path.exists(jobsDir):
    os.makedirs(jobsDir)


if args.genome=='mm9':
	args.genome =' /lustre/beagle2/ReferenceSequences/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome '
elif args.genome=='mm10':
	args.genome =' /lustre/beagle2/ReferenceSequences/Mus_musculus/GENCODE/mm10/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome '
elif args.genome=='hg19':
	args.genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/GRCh37.p13.genome '
elif args.genome == 'hg38':
	args.genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh38.primary_assembly.genome '
else:
	print("specify mm9, mm10, hg19, or hg38")
	
	
FilesList= glob.glob(os.path.join(args.fastqDir, '*.fastq'))
for fs in FilesList:
    f=fs.strip().split('/')[-1].strip().split('.')[0]
    fout=open(jobsDir+f+'.pbs','w')
    fout.write('#!/bin/bash' + ' \n')
    fout.write('#PBS -N toph' + '\n')
    fout.write('#PBS -S /bin/bash' + '\n')
    fout.write('#PBS -l mppwidth='+args.processors+ '\n')
    fout.write('#PBS -l walltime='+args.walltime+':00:00' + '\n')
    fout.write('#PBS -j oe  \n')
    fout.write('#PBS -o ' + jobsDir +  '/'+f+'.log'+'\n')
    fout.write('#PBS -e ' + jobsDir  + '/'+f+'.err' + '\n')
    fout.write('#PBS -l mppwidth=32' + '\n')
    fout.write('cd $PBS_O_WORKDIR'+'\n')
    fout.write('aprun -n 1 -N 1 -d ' +args.processors+ '  /lustre/beagle2/cperez5/soft/shellScripts/bowtie2_runner.sh  -o ' +args.out+ ' -g '+args.genome+ '  -f ' +fs+  ' -c ' +args.processors+ '  -s '+f+  ' -m ' +args.mapping+ ' \n')
    fout.close()
    time.sleep(1)
    os.system('qsub '+jobsDir+ f + '.pbs')












	
