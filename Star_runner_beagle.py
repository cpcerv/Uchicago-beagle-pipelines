#!/bin/bash

##### WRITE PBS SCRIPT AND CALL STAR #########
#	This python script will write a .pbs for every fastq file      #
#	in a directory.  It calls a master tophat-runner shell script #
#	for aprun.  n and N for aprun will always be 1.                 #
#	Open source. By:Carlos PerezCervantes. Moskowitz lab # 
#####################################################

import os, datetime, glob, time, argparse


###########################################################
##
##  define parser and parameters for user input
###########################################################

parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string) # allows a directory to be used as input
	return os.path.abspath(string)


parser.add_argument("-d", "--fastqDir", type=readable_dir,
		help="directory of Fastq files")
				  
parser.add_argument("-o", "--out", type=str,
		help="specify output destination")
				  
parser.add_argument("-g", "--genome", type=str,default='mm10',
		help="specify genome: mm10, hg19, hg38")
				  
parser.add_argument("-p", "--processors", type=str, default='32',
		help="specify number of processors for computation node (a mutiple of 32) and tophat")              

parser.add_argument("-w", "--walltime", type=str, default='02', 
		help="specify running time in hours, i.e 02, 4, 16")

args = parser.parse_args()
args.out = os.path.abspath(args.out)

##########################################3
##
##   create log and pbs directory
##
###########################################

if not os.path.exists(args.out):
    os.makedirs(args.out)

jobsDir = args.out+'/log/'

if not os.path.exists(jobsDir):
    os.makedirs(jobsDir)

###########################################
##
##  Define genome locations from input
###########################################


if args.genome=='mm10':
	args.genome ='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Sequence/STARgenome/ '
elif args.genome=='hg19':
	args.genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/STARgenome_hg19/ '
elif args.genome == 'hg38':
	args.genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/STARgenome/  '
else:
	print("specify mm10, hg19, or hg38")
	

##########################################################
##
##  create PBS and define job dependencies
##
##########################################################

	
FilesList= glob.glob(os.path.join(args.fastqDir, '*.fastq.gz'))
for fs in FilesList:
    f=fs.strip().split('/')[-1].strip().split('.')[0]
    if not os.path.exists(args.out+'/'+f):
        os.makedirs(args.out+'/'+f)
    fout=open(jobsDir+f+'.pbs','w')
    fout.write('#!/bin/bash' + ' \n')
    fout.write('#PBS -N s'+f.replace("-","")+ '  ' + '\n')
    fout.write('#PBS -S /bin/bash' + '\n')
    fout.write('#PBS -l mppwidth='+args.processors+ '\n')
    fout.write('#PBS -l walltime='+args.walltime+':00:00' + '\n')
    fout.write('#PBS -j oe  \n')
    fout.write('#PBS -o ' + jobsDir +  '/'+f+'.log'+'\n')
    fout.write('#PBS -e ' + jobsDir  + '/'+f+'.err' + '\n')
    fout.write('#PBS -l mppwidth=32' + '\n')
    fout.write('cd $PBS_O_WORKDIR'+'\n')
    fout.write('if [ ! $(module list -t 2>&1 | grep PrgEnv-gnu) ]; then '+ ' \n')
    fout.write('module swap PrgEnv-cray PrgEnv-gnu'+'\n')
    fout.write('fi' + '\n')
    fout.write('aprun -n 1 -N 1 -d ' +args.processors )
    fout.write(" /soft/STAR/gnu/bin/STAR --runThreadN "+args.processors+" --genomeDir "+args.genome+" --readFilesIn "+fs+" --readFilesCommand zcat --outSAMstrandField  --outFileNamePrefix "+args.out+"/"+f+"/  --outSAMtype BAM SortedByCoordinate")
    fout.close()
    time.sleep(1)
    #os.system('qsub '+jobsDir+ f + '.pbs')




