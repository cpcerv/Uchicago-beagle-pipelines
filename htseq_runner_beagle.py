#!/bin/bash
import os, datetime, glob, time, argparse

################################################################
#  	open source HTseq runner version 0.1 for beagle torque cluster         #
#    			Carlos Perez-Cervantes; Moskowitz Lab			      #
# 	Input a directory of bam or sam files and a gtf and output count files   #
#													     #
###############################################################


#####################################################\
# 1. define functions # \
#####################################################\

parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string)
	return os.path.abspath(string)

#####################################################\
# 2. create argument list  and check conditions# \
#####################################################\


parser.add_argument("-d", "--bamDir", type=readable_dir,
		help="directory of Bam files")#, default=os.getcwd()+'/FastQ')
				  
parser.add_argument("-o", "--out", type=str,
		help="specify output destination")

parser.add_argument("-i", "--id", type=str,default='gene_id',
		help="specify name of gene id column to use, e.g gene_name, gene_id")
				  
parser.add_argument("-g", "--gtf", type=str, default=' mm10',
		help="specify gtf to use : mm9,  mm10, hg19, hg38 or custom  ")
				  
parser.add_argument("-p", "--processors", type=str, default='32',
		help="specify number of processors, a multiple of 32")              

parser.add_argument("-w", "--walltime", type=str, default='03', 
		help="specify running time in hours, i.e 02, 04, 16")

parser.add_argument("-s", "--stranded", type=str, default='reverse', 
		help="specify strandedness, options are 'yes','no', or 'reverse' ")

parser.add_argument("-m", "--mode", type=str, default='union', 
		help="specify HTseq counting mode ie  'union', 'intersection-strict', 'intersection-nonempty' ")
		
args = parser.parse_args()
args.out = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)

jobsDir = args.out+'/log/'

if not os.path.exists(jobsDir):
    os.makedirs(jobsDir)

if args.gtf=='mm9':

	args.gtf =' /lustre/beagle2/ReferenceSequences/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome.fa '
	gtf= '/lustre/beagle2/ReferenceSequences/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-current/Genes/genes.gtf'


elif args.gtf=='mm10':

	args.gtf =' /lustre/beagle2/ReferenceSequences/Mus_musculus/GENCODE/mm10/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome.fa '
	gtf='/lustre/beagle2/ReferenceSequences/Mus_musculus/GENCODE/mm10/release_M9/gencode.vM9.annotation.gtf'
	bed='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Annotation/gencode.vM9.gene_protein_coding.gtf'

elif args.gtf=='hg19':

	args.gtf =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/GRCh37.p13.genome.fa '
	gtf='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v19.annotation.gtf'
	bed='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v19.gene_protein_coding.gtf'

elif args.gtf == 'hg38':

	args.gtf =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh38.primary_assembly.genome.fa '
	gtf='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v25.primary_assembly.annotation.gtf'
	bed='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v25.gene_protein_coding.gtf'

elif args.gtf == 'custom':
	args.gtf = raw_input("Enter gtf i.e /path/merged.gtf"+'\n')
	args.gtf = os.path.abspath(args.gtf)
else:
	print("specify mm9, mm10, mm10, hg19, hg38 or custom")



#####################################################\
# 3. Print Resource Manager Directives to .pbs file # 
#####################################################\

FilesList = glob.glob(os.path.join(args.bamDir, '*.bam'))
for fs in FilesList:
	f=fs.strip().split('/')[-1].strip().split('_')[0]
	if fs.endswith(".bam"):
		filetype = 'bam'
	elif fs.endswith(".sam"):
		filetype == ".sam"
	else: 
		print("No bam or sam files here")
	fout = open(jobsDir + f + '.pbs','w')
	fout.write('#PBS -S /bin/bash' + '\n')
	fout.write('#PBS -N htseq_'+f.replace("-","") + ' \n')
	fout.write('#PBS -l mppwidth='+args.processors+ '\n')
	fout.write('#PBS -l walltime='+args.walltime+':00:00' + '\n')
	fout.write('#PBS -o '+jobsDir+f+'.log ' +  '\n')
	fout.write('#PBS -e '+jobsDir+f+'.err ' +   '\n')
	fout.write('' + '\n')
	fout.write('. /opt/modules/default/init/bash'+'\n')
	fout.write('if [ ! $(module list -t 2>&1 | grep PrgEnv-gnu) ]; then'+'\n')
	fout.write('module swap PrgEnv-cray PrgEnv-gnu'+'\n')
	fout.write('module load python/2.7.6-vanilla ' + '\n')
	fout.write('fi'+'\n')
	fout.write('cd $PBS_O_WORKDIR'+'\n')
	fout.write(''+'\n')
	fout.write('aprun -n 1 -N 1 -d 32 ' )
#####################################\
# 4. Print commands to .pbs file # 
#####################################\
	fout.write(' htseq-count -f ' +filetype+  ' -s ' +args.stranded+  ' -r pos -i ' + args.id+ ' ' + fs + ' ' + args.gtf+ ' -m ' +args.mode)
	fout.write(' > ' + args.out + '/' +f + '.counts')
	fout.close()
	os.system('qsub '+jobsDir+f+ '.pbs ')
#	os.system(f.strip().split('-')[1]+'=$(qsub -W depend=afterok:$4435972.sdb '+jobsDir+ f + '.pbs)')
#	os.system('echo ' +f.strip().split('-')[1])    













