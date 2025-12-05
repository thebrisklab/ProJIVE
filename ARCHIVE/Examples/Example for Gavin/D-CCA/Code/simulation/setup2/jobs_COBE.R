#dqshtc

flist=  list.files("/scratch/hshu/paper1_revision/setup2/revision2", pattern = "*_p900_angle45_nstd1_COBE.txt")
fdone=as.integer(sub('.*siml2_seed*(.*?) *_p900_angle45_nstd1_COBE.txt.*','\\1',flist))
#for(j in 1:1000)
for(j in setdiff(1:1000,fdone))
{
  script=NULL
  script[1]=sprintf("#PBS -N COBE2")
  script[2] ="#PBS -l procs=1,pmem=2gb\n
  #PBS -l walltime=2:00:00"      
  
  script[3] = "cd /scratch/hshu/paper1_revision/setup2/revision2\n
  module load matlab/2014b\n"
  
  script[4] = sprintf("matlab -nodisplay -nodesktop -singleCompThread -r 'simulation2_COBE(%d)'",j)
  
  write.table(script,file = sprintf("/scratch/hshu/paper1_revision/setup2/jobs/job_simulation2_COBE%d.pbs",j),row.names=F,col.names = F,quote=F)
  
  system(sprintf("qsub /scratch/hshu/paper1_revision/setup2/jobs/job_simulation2_COBE%d.pbs",j),intern=F)
}    




#eagle
flist=  list.files("/rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2", pattern = "*_p900_angle45_nstd1_COBE.txt")
fdone=as.integer(sub('.*siml2_seed*(.*?) *_p900_angle45_nstd1_COBE.txt.*','\\1',flist))

    for(j in 1:1000)
    {
      script=NULL
      script[1]=sprintf("#BSUB -J COBE2")
      script[2] ="#BSUB -W 2:00\n
#BSUB -q medium\n
#BSUB -n 1\n
#BSUB -M 8192\n
#BSUB -R rusage[mem=8192]\n
#BSUB -N \n"      
      
      script[3] = "#BSUB -o /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2/data/output\n
#BSUB -e /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2/data/output\n
#BSUB -cwd /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2\n
      module load Matlab\n"
      
      script[4] = sprintf("matlab -nodisplay -nodesktop -singleCompThread -r 'simulation2_COBE(%d)'",j)
      
      write.table(script,file = sprintf("/rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2/data/jobs/job_simulation1_COBE%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("bsub < /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2/data/jobs/job_simulation1_COBE%d.pbs",j),intern=F)
    }
