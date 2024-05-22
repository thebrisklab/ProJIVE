#python3/3.5.1-0-anaconda
    for(j in 1:1000)
    {
      script=NULL
      script[1]=sprintf("#BSUB -J AJIVE2")
      script[2] ="#BSUB -W 2:00\n
#BSUB -q medium\n
#BSUB -n 1\n
#BSUB -M 8192\n
#BSUB -R rusage[mem=8192]\n
#BSUB -N \n"      
      
      script[3] = "#BSUB -o /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data/output\n
#BSUB -e /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data/output\n
#BSUB -cwd /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data\n
module load Matlab\n"
      
      script[4] = sprintf("matlab -nodisplay -nodesktop -singleCompThread -r 'simulation2_AJIVE(%d)'",j)
      
      write.table(script,file = sprintf("/rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data/jobs/job_simulation1_AJIVE%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("bsub < /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data/jobs/job_simulation1_AJIVE%d.pbs",j),intern=F)
    }

    
    
#dqshtc
    
    flist=  list.files("/scratch/hshu/paper1_revision/setup2/revision2", pattern = "*_p900_angle45_nstd1_AJIVE.txt")
    fdone=as.integer(sub('.*siml2_seed*(.*?) *_p900_angle45_nstd1_AJIVE.txt.*','\\1',flist))
    #for(j in 1:1000)
    for(j in setdiff(1:1000,fdone))
    {
      script=NULL
      script[1]=sprintf("#PBS -N AJIVE2")
      script[2] ="#PBS -l procs=1,pmem=2gb\n
#PBS -l walltime=2:00:00"      
      
      script[3] = "cd /scratch/hshu/paper1_revision/setup2/revision2\n
      module load matlab/2014b\n"
      
      script[4] = sprintf("matlab -nodisplay -nodesktop -singleCompThread -r 'simulation2_AJIVE(%d)'",j)
      
      write.table(script,file = sprintf("/scratch/hshu/paper1_revision/setup2/jobs/job_simulation2_AJIVE%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("qsub /scratch/hshu/paper1_revision/setup2/jobs/job_simulation2_AJIVE%d.pbs",j),intern=F)
    }    
    
    
    