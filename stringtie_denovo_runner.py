import os, datetime, glob, time, argparse, shutil
parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string)
	return os.path.abspath(string)

parser.add_argument("-d", "--bamDir", type=readable_dir,
		help="directory of bam files")


parser.add_argument("-o", "--out", type=str,
		help="specify output destination i.e ../stringtie-out will create folder if needed")


parser.add_argument("-g", "--genome", type=str,default='mm10',
		help="specify genome: mm9, mm10, hg19, hg38")

				  
parser.add_argument("-p", "--processors", type=str, default='32',
		help="specify number of processors for computation node (a mutiple of 32) and cufflinks")              


parser.add_argument("-w", "--walltime", type=str, default='16', 
		help="specify running time in hours, i.e 02, 04, 16")

parser.add_argument("-s", "--strand", type=str, default='rf', 
        help="specify strandedness, rf for  fr-firststrand , fr for fr-secondstrand")

		
args = parser.parse_args()
args.out = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)

jobsDir = args.out+'/log/'

if not os.path.exists(jobsDir):
    os.makedirs(jobsDir)

if args.genome=='mm9':

	args.genome =' /lustre/beagle2/ReferenceSequences/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome.fa '
	gtf= '/lustre/beagle2/ReferenceSequences/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-current/Genes/genes.gtf'


elif args.genome=='mm10':

	args.genome =' /lustre/beagle2/ReferenceSequences/Mus_musculus/GENCODE/mm10/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome.fa '
	gtf='/lustre/beagle2/ReferenceSequences/Mus_musculus/GENCODE/mm10/release_M9/gencode.vM9.annotation.gtf'


elif args.genome=='hg19':

	args.genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/GRCh37.p13.genome.fa '
	gtf='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v19.annotation.gtf'


elif args.genome == 'hg38':

	args.genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh38.primary_assembly.genome.fa '
	gtf='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v25.primary_assembly.annotation.gtf'


else:
	print("specify mm9, mm10, hg19, or hg38")
	
	
FilesList = glob.glob(os.path.join(args.bamDir,  '*[!.bai].bam'))
for fs in FilesList:
    f=fs.strip().split('/')[-1].strip().split('.')[0]
    if not os.path.exists(args.out+'/'+f):
        os.makedirs(args.out+'/'+f)
    fout=open(jobsDir+f+'.pbs','w')
    fout.write('#!/bin/bash' + ' \n')
    fout.write('#PBS -N stie'+f.replace("-","") + '\n')
    fout.write('#PBS -S /bin/bash' + '\n')
    fout.write('#PBS -l mppwidth='+args.processors+ '\n')
    fout.write('#PBS -l walltime='+args.walltime+':00:00' + '\n')
    fout.write('#PBS -j oe  \n')
    fout.write('#PBS -o ' + jobsDir +  '/'+f+'.log'+'\n')
    fout.write('#PBS -e ' + jobsDir  + '/'+f+'.err' + '\n')
    fout.write('cd $PBS_O_WORKDIR'+'\n')
    fout.write('if [ ! $(module list -t 2>&1 | grep PrgEnv-gnu) ]; then '+ ' \n')
    fout.write('module swap PrgEnv-cray PrgEnv-gnu'+'\n')
    fout.write('moduleload stringtie'+'\n')
    fout.write('aprun -n 1 -N 1 -d ' +args.processors)
    fout.write(' stringtie ' +fs+ ' -l '+f+ ' --' +args.strand)
    fout.write(' -m 50 -p '+args.processors)
    fout.write(' -o ' +args.out+'/'+f+'/transcripts.gtf')
    fout.write(' '+'\n')
    fout.close()
    time.sleep(1)
    os.system('qsub '+jobsDir+ f + '.pbs')




