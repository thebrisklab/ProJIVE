#flist=  list.files("/scratch/hshu/paper1_revision/setup1/revision2", pattern = "*_p900_angle45_nstd1_OnPLS.txt")
#fdone=as.integer(sub('.*siml1_seed*(.*?) *_p900_angle45_nstd1_OnPLS.txt.*','\\1',flist))

    #dqshtc #python3/3.5.2-0-anaconda
    #for(j in setdiff(1:1000,fdone))
    for(j in 1:1000)
    {
      script=NULL
      script[1]=sprintf("#PBS -N OnPLS1")
      script[2] ="#PBS -l procs=1,pmem=2gb\n
#PBS -l walltime=120:00:00\n"      
      
      script[3] = "cd /scratch/hshu/paper1_revision/setup1/revision2\n
      setenv PATH /workspace/hshu/anaconda3/bin:$PATH\n
      setenv PYTHONPATH /workspace/hshu/anaconda3/bin:$PATH\n
      setenv MKL_NUM_THREADS 1\n"
      
      script[4] = sprintf("python simulation1_OnPLS.py %d",j)
      
      write.table(script,file = sprintf("/scratch/hshu/paper1_revision/setup1/jobs/job_simulation1_OnPLS%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("qsub /scratch/hshu/paper1_revision/setup1/jobs/job_simulation1_OnPLS%d.pbs",j),intern=F)
    }    
    
    
