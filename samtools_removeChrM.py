import os, datetime, glob, time, argparse, shutil
parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string)
	return os.path.abspath(string)

parser.add_argument("-d", "--bamDir", type=readable_dir,
		help="directory of bam files")

parser.add_argument("-o", "--out", type=str,
		help="specify output destination")

args = parser.parse_args()

args.out = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)

jobsDir = args.out+'/log/'

if not os.path.exists(jobsDir):
    os.makedirs(jobsDir)


FilesList= glob.glob(os.path.join(args.bamDir, '*[!.bai].bam'))

for fs in FilesList:
    f=fs.strip().split('/')[-1].strip().split('_')[0]
    fout=open(jobsDir+f+'.pbs','w')
    fout.write('#!/bin/bash'+ ' \n')
    fout.write(''+'\n')
    fout.write('#PBS -S /bin/bash' + '\n')
    fout.write('#PBS -N ' + f + '\n')
    fout.write('#PBS -l mppwidth=32 \n')
    fout.write('#PBS -l walltime=10:00:00 ' + '\n')
    fout.write('. /opt/modules/default/init/bash'+'\n')
    fout.write('if [ ! $(module list -t 2>&1 | grep PrgEnv-gnu) ]; then'+'\n')
    fout.write('module swap PrgEnv-cray PrgEnv-gnu'+'\n')
    fout.write('module load samtools/1.5 ' + '\n')
    fout.write('fi'+'\n')
    fout.write('cd $PBS_O_WORKDIR'+'\n')
    fout.write(''+'\n')
    fout.write('aprun -n 1 -N 1 -d 32' )
    fout.write('samtools idxstats ' + fs + ' ' + '| cut -f 1 | grep -v chrM |' + ' ' + 'xargs samtools view -q 30  -b ' + fs + ' > ' +  args.out +'/' + f  + ' \n')
    fout.close()
    time.sleep(1)
    os.system('qsub '+jobsDir+ f + '.pbs')
