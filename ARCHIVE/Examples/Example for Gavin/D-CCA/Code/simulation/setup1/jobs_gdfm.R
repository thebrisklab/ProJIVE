#eagle
    for(j in 1:1000)
    {
      script=NULL
      script[1]=sprintf("#BSUB -J gdfm1")
      script[2] ="#BSUB -W 2:00\n
#BSUB -q medium\n
#BSUB -n 1\n
#BSUB -M 8192\n
#BSUB -R rusage[mem=8192]\n
#BSUB -N \n"      
      
      script[3] = "#BSUB -o /rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/data/output\n
#BSUB -e /rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/data/output\n
#BSUB -cwd /rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/data\n
      module load Matlab\n"
      
      script[4] = sprintf("matlab -nodisplay -nodesktop -singleCompThread -r 'simulation1_gdfm(%d)'",j)
      
      write.table(script,file = sprintf("/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/data/jobs/job_simulation1_gdfm%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("bsub < /rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/data/jobs/job_simulation1_gdfm%d.pbs",j),intern=F)
    }


    
#dqshtc
    flist=  list.files("/scratch/hshu/paper1_revision/setup1/revision2", pattern = "*_p900_angle45_nstd1_gdfm_nlogP_2.txt")
    fdone=as.integer(sub('.*siml1_seed*(.*?) *_p900_angle45_nstd1_gdfm_nlogP_2.txt.*','\\1',flist))#933
  #  for(j in 1:25)
    for(j in setdiff(1:1000,fdone))
    {
      script=NULL
      script[1]=sprintf("#PBS -N gdfm1")
      script[2] ="#PBS -l procs=1,pmem=2gb\n
#PBS -l walltime=1:00:00"      
      
      script[3] = "cd /scratch/hshu/paper1_revision/setup1/revision2\n
      module load matlab/2014b\n"
      
      script[4] = sprintf("matlab -nodisplay -nodesktop -singleCompThread -r 'simulation1_gdfm(%d)'",j)
      
      write.table(script,file = sprintf("/scratch/hshu/paper1_revision/setup1/jobs/job_simulation1_gdfm%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("qsub /scratch/hshu/paper1_revision/setup1/jobs/job_simulation1_gdfm%d.pbs",j),intern=F)
    }    
    
    
    
    
    
    
    