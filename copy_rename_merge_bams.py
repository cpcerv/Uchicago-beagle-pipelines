#!/bin/bash

import os, datetime, glob, time, argparse

parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string)
	return os.path.abspath(string)


parser.add_argument("-d", "--tophatDir", type=readable_dir,
		help="directory of tophat files")
				  
parser.add_argument("-o", "--out", type=str,
		help="specify output destination i.e '-o output' ")
				  			
parser.add_argument("-p", "--processors", type=str, default='32',
		help="specify number of processors for computation node (a mutiple of 32) and tophat")              

parser.add_argument("-w", "--walltime", type=str, default='8', 
		help="specify running time in hours, i.e 02, 4, 16")


parser.add_argument("-s", "--submit", type=bool, default=True, 
        help="submit PBS True or False")


args = parser.parse_args()
args.out = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)

jobsDir = args.out+'/log/'

if not os.path.exists(jobsDir):
    os.makedirs(jobsDir)

FilesList= glob.glob(os.path.join(args.tophatDir,'*[!log]'))

samples = {}
for fastqzfile in FilesList:
    f=fastqzfile.strip().split('/')[-1].strip().split('.')[0]
    name= f.strip().split('_')[0]
    if not ("Undetermined" in f):
        if not (name in samples):
            samples[name]=[]
        src = args.tophatDir +'/'+ f +'/accepted_hits.bam'
        samples[name].append(src)

        
for name in samples:
    dst = args.out +'/'+ name +'.bam'
    if len(samples[name]) == 1:
        os.system('cp ' + samples[name][0] + ' ' + dst)
    if len(samples[name]) > 1:
        file_list = ' '.join(samples[name])
        fout=open(jobsDir+name+'.pbs','w')
        fout.write('#!/bin/bash' + ' \n')
        fout.write('#PBS -S /bin/bash' + '\n')
        fout.write('#PBS -N '+name.replace("-","")+ '\n')
        fout.write('#PBS -l mppwidth='+args.processors+ '\n')
        fout.write('#PBS -l walltime='+args.walltime+':00:00' + '\n')
        fout.write('#PBS -j oe  \n')
        fout.write('#PBS -o ' + jobsDir + '/'+ name + '.log'+'\n')
        fout.write('#PBS -e ' + jobsDir + '/'+ name + '.err' + '\n')
        fout.write('. /opt/modules/default/init/bash'+'\n')
        fout.write('if [ ! $(module list -t 2>&1 | grep PrgEnv-gnu) ]; then'+'\n')
        fout.write('module swap PrgEnv-cray PrgEnv-gnu'+'\n')
        fout.write('module load samtools/1.5 ' + '\n')
        fout.write('fi'+'\n')
        fout.write('cd $PBS_O_WORKDIR'+'\n')
        fout.write(''+'\n')
        fout.write('aprun -n 1 -N 1 -d ' +args.processors)
        fout.write(' samtools merge -f ' + dst + ' ' + file_list + '\n')
        fout.close()
    
        time.sleep(1)
        if args.submit==True:
            os.system('qsub '+jobsDir+ name + '.pbs')
        else:
            break
        

