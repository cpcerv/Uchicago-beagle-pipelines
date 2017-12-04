#!/usr/bin/env python
import os, datetime, glob, time, argparse, subprocess





#################### RNAaseq pipeline for PolyA++ reads ################################
##      This python script will write a .pbs for every fastq file                     #
##      in a directory.  it will call tophat,samtools and HTseq.                      #
##      These will be output to the same directory branch as FastQ                    #
##      assumes single-end reads                                                      #
##      n and N for aprun will always be 1                                            #
##      Open source. By:Carlos Perez-Cervantes. Moskowitz lab                         # 
########################################################################################




parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string)
	return os.path.abspath(string)


parser.add_argument("-d", "--fastqDir", type=readable_dir,
		help="directory of Fastq files")
				  
parser.add_argument("-o", "--out", type=str,
		help="specify output destination for all steps in the pipeline")
				  
parser.add_argument("-g", "--genome", type=str,default='mm10',
		help="specify genome: mm9, mm10, hg19, hg38")
				  
parser.add_argument("-p", "--processors", type=str, default='32',
		help="specify number of processors for computation node (a mutiple of 32) and tophat") 

parser.add_argument("-w", "--walltime", type=str, default='8', 
		help="specify running time in hours, i.e 02, 4, 16")

parser.add_argument("-l", "--libraryType", type=str, default='firststrand', 
		help="specify library type, options are 'firststrand','secondstrand', or 'unstranded' ")

parser.add_argument("-a", "--aligner", type=str, default='tophat',
    help="specify which aligner to use.  options: tophat, or STAR")
		
args = parser.parse_args()
parent = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)


###############################################################################
##
##  assign annotation from user input
##
###############################################################################


if args.genome=='mm9':
    genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/mm9/Bowtie2Index/genome '
    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/mm9/Bowtie2Index/genome.fa '
    gtf= '/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/mm9/genes.gtf '

elif args.genome=='mm10':

    
    if args.aligner ==  "STAR":
      genome_index = '/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Sequence/STARgenome/ '
    else:
      genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome '
      genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome.fa '
    gtf=' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Annotation/gencode.vM9.annotation.gtf'


elif args.genome=='hg19':

    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/GRCh37.p13.genome.fa '
    gtf='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v19.annotation.gtf'
    genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/GRCh37.p13.genome '


elif args.genome == 'hg38':
    if args.aligner ==  "STAR":
      genome_index = '/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/STARgenome/ '
    else:
      genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh38.primary_assembly.genome.fa '
      genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh38.primary_assembly.genome '
    gtf='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v25.primary_assembly.annotation.gtf'

else:
    print("specify mm9, mm10, hg19, or hg38")


##############################################################################
##  make directories for output
##
##############################################################################

parentList= glob.glob(os.path.join(args.fastqDir, '[!Undetermined]*.fastq.gz'))


mainseed = list(set([x.strip().split('/')[-1].strip().split('_')[0] for x in parentList]))





if args.aligner == "tophat":
  alignDir = parent+'/tophat_'+args.genome+'/'
  if not os.path.exists(alignDir):
    os.makedirs(alignDir)


elif args.aligner == "STAR":
  alignDir = parent+'/STAR_'+args.genome+'/'
  if not os.path.exists(alignDir):
    os.makedirs(alignDir)


alignDir = os.path.abspath(alignDir)


#all_bam = parent+'/all_bam_'+args.genome+'/'


#if not os.path.exists(all_bam):
#   	os.makedirs(all_bam)



#all_bam = os.path.abspath(all_bam)



cleaned_bam = parent+'/cleaned_bam_'+args.genome+'/'

if not os.path.exists(cleaned_bam):
   	os.makedirs(cleaned_bam)

cleaned_bam = os.path.abspath(cleaned_bam)



htseq = parent+'/htseq_'+args.genome+'/'

if not os.path.exists(htseq):
   	os.makedirs(htseq)

htseq = os.path.abspath(htseq)



##############################################################################
##
##   make pbs and log directories
##
##############################################################################


if not os.path.exists(alignDir+'/log/'):
    os.makedirs(alignDir+'/log/')

algPBS = alignDir+'/log/'

if not os.path.exists(cleaned_bam+'/log/'):
        os.makedirs(cleaned_bam+'/log/')

qcjob = cleaned_bam+'/log/'


if not os.path.exists(htseq+'/log/'):
        os.makedirs(htseq+'/log/')

htPBS = htseq+'/log/'



##############################################################################
##  define function that will be looped to write many 
##  PBS files
##############################################################################


def pbs_writer(pbsDir, seed, processors, walltime, shell):
    pbsf=open(pbsDir+seed+'.pbs','w')
    pbsf.write('#!/bin/bash' + ' \n')
    pbsf.write('#PBS -N ' +seed+ '\n')
    pbsf.write('#PBS -S /bin/bash' +'\n')
    pbsf.write('#PBS -l mppwidth='+processors+ '\n')
    pbsf.write('#PBS -l walltime='+walltime+':00:00' + '\n')
    pbsf.write('#PBS -j oe  \n')
    pbsf.write('#PBS -o ' + pbsDir +seed+'.log'+'\n')
    pbsf.write('#PBS -e ' + pbsDir +seed+'.err' + '\n')
    pbsf.write('./opt/modules/default/init/bash '+'\n')
    pbsf.write('cd $PBS_O_WORKDIR'+'\n')
    pbsf.write('module swap PrgEnv-cray PrgEnv-gnu'+'\n')
    pbsf.write('module load samtools/1.5 '+'\n')
    pbsf.write('aprun -n 1 -N 1 -d ' +processors+' ' +shell)
    pbsf.close()


##############################################################################
##
##  run aligner on FastQ file directory. ##
##############################################################################
 
if args.aligner == "STAR":
  for fs in parentList:
    f=fs.strip().split('/')[-1].strip().split('.')[0]
    if not os.path.exists(alignDir+'/'+f.strip().split('_')[0]):
      os.makedirs(alignDir+'/'+f.strip().split('_')[0])
    pbs_writer(pbsDir = algPBS, seed=f.strip().split('_')[0],
    processors=str(args.processors), walltime=str(args.walltime),
    shell = '/lustre/beagle2/cperez5/soft/shellScripts/STAR_runner.sh  -o ' +alignDir+ '  -f ' +fs+  ' -p ' +str(args.processors)+ '  -s '+f.strip().split('_')[0]+' -c zcat -g ' +genome_index) 
  

else:
  for fs in parentList:
    f=fs.strip().split('/')[-1].strip().split('.')[0]
    pbs_writer(pbsDir = algPBS, seed=f.replace("-",""),
    processors=str(args.processors), walltime=str(args.walltime),
    shell = '/lustre/beagle2/cperez5/soft/shellScripts/tophat2_runner.sh  -o ' +alignDir+ ' -g '+genome_index+ '  -f ' +fs+  ' -c ' +str(args.processors)+ '  -s '+f)




##############################################################################
##	Write pbs for Index and clean bam files
##
##############################################################################


for fs in mainseed:
    f=fs.strip().split('/')[-1].strip().split('_')[0]
    if args.aligner == "STAR":
      pbs_writer(pbsDir = qcjob,  seed='qc_'+f, walltime = str(args.walltime),
      processors=str(args.processors), shell = '/lustre/beagle2/cperez5/soft/shellScripts/samtoolsQC_runner.sh -o ' +cleaned_bam+ ' -f ' +alignDir+'/'+f+ '/Aligned.sortedByCoord.out.bam -s '+f)
    else:
      pbs_writer(pbsDir = qcjob,  seed='qc_'+f, walltime = str(args.walltime),
      processors=str(args.processors), shell = '/lustre/beagle2/cperez5/soft/shellScripts/samtoolsQC_runner.sh -o ' +cleaned_bam+ ' -f ' +alignDir+'/'+f+ '/accepted_hits.bam -s '+f)



##############################################################################
##	Write HTseq pbs scripts
##
##############################################################################


for fs in mainseed:
    pbs_writer(pbsDir= htPBS, seed='ht_'+fs, walltime=str(args.walltime),
    processors=str(args.processors),
    shell = '/lustre/beagle2/cperez5/soft/shellScripts/htseq_runner.sh -t bam -i gene_name -r reverse -a ' +gtf + ' -o ' +htseq+ ' -s ' +fs+ ' -f '+cleaned_bam+'/'+fs+'.bam')



######################################
##
##  Get mapping stats and plot with R
##
######################################


#pbs_writer(pbsDir = algPBS, seed = 'tophSummry' , walltime= '01', processors= str(args.processors), shell=' python /lustre/beagle2/cperez5/soft/python2Scripts/get_tophat_summary.py -d ' +alignDir)

pbs_writer(pbsDir = qcjob,  seed='bamsummary', walltime ='01', processors=str(args.processors), shell = 'python /lustre/beagle2/cperez5/soft/python2Scripts/samtools_bamsummary_local.py -d ' +cleaned_bam)

pbs_writer(pbsDir = htPBS, seed = 'htsqSummry' , walltime= '01', processors= str(args.processors), shell='python /lustre/beagle2/cperez5/soft/python2Scripts/htseq_summary.py -d '+htseq)

pbs_writer(pbsDir = htPBS, seed = 'plotRline' , walltime= '01', processors= str(args.processors), shell='R CMD BATCH /lustre/beagle2/cperez5/soft/R_scripts/plot_line.py '+algPBS+'tophSummry.pbs' + qcjob+' bamsummary.pbs '+htPBS+'htsqSummry.pbs ' )



##############################################################################
##
##  Run all qsub commands
##
##############################################################################




htslist=[]

if os.path.isfile(args.out+"/submitALL.sh"):
  os.remove(args.out+"/submitALL.sh")
with open(args.out+'/submitALL.sh', "a") as runnit:
    runnit.write("#!/bin/bash\n \n \n")


for name in mainseed:
  alcmd = 'al'+name.replace("-","")+'=$(qsub '+algPBS+name+ '.pbs)\necho -e $al'+name.replace("-","")+ ' \n'
  qccm='qc'+name.replace("-","")+'=$(qsub -W depend=afterok:$al'+name.replace("-","")+ ' '+qcjob+'qc_'+name+'.pbs) \n echo -e $qc'+name.replace("-","")+ ' \n'
  summ2='bamsummary=$(qsub -W depend=afterok:+ $qc'+name.replace("-","") +' '+'bamsummary.pbs) \n echo -e $bamsummary.pbs\n'
  htscmd='hts'+name.replace("-","")+'=$(qsub -W depend=afterok:$qc'+name.replace("-","") +' '+htPBS+'ht_'+name+ '.pbs) \n echo -e $hts'+name.replace("-","") + ' \n'
  #summ3='htsum'='$(qsub -W depend=afterok:'+":".join([str(i) for i in htslist]) +' htsqSummry.pbs) \n echo -e $htsum\n'
  with open(args.out+'/submitALL.sh', "a") as runnit:
    runnit.write(alcmd+ '\nsleep 2\n\n\n\n')
    runnit.write(qccm+ '\nsleep 2\n\n\n\n')
    runnit.write(htscmd+ '\nsleep 2\n\n\n\n')
    #runnit.write(summ1+ '\nsleep 2\n\n\n\n')
    #runnit.write(summ2+ '\nsleep 2\n\n\n\n')
    #runnit.write(summ3+ '\nsleep 2\n\n\n\n')
    runnit.close()

os.system('chmod u+x ' +args.out+"/submitALL.sh")
os.system(args.out+"/submitALL.sh")

