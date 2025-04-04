#python3/3.5.1-0-anaconda
    for(j in 1:1000)
    {
      script=NULL
      script[1]=sprintf("#BSUB -J OnPLS2")
      script[2] ="#BSUB -W 40:00\n
#BSUB -q long\n
#BSUB -n 1\n
#BSUB -M 8192\n
#BSUB -R rusage[mem=8192]\n
#BSUB -N \n"      
      
      script[3] = "#BSUB -o /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data/output\n
#BSUB -e /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data/output\n
#BSUB -cwd /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data\n"
      
      script[4] = sprintf("python simulation2_OnPLS.py %d",j)
      
      write.table(script,file = sprintf("/rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data/jobs/job_simulation1_OnPLS%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("bsub < /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision1/data/jobs/job_simulation1_OnPLS%d.pbs",j),intern=F)
    }
    
    
    #dqshtc #python3/3.5.2-0-anaconda
    #flist=  list.files("/scratch/hshu/paper1_revision/setup2/revision2", pattern = "*_p900_angle45_nstd1_OnPLS.txt")
    #fdone=as.integer(sub('.*siml2_seed*(.*?) *_p900_angle45_nstd1_OnPLS.txt.*','\\1',flist))
    
    for(j in 1:1000)
    {
      script=NULL
      script[1]=sprintf("#PBS -N OnPLS2")
      script[2] ="#PBS -l procs=1,pmem=2gb\n
#PBS -l walltime=40:00:00\n"      
      
      script[3] = "cd /scratch/hshu/paper1_revision/setup2/revision2\n
      setenv PATH /workspace/hshu/anaconda3/bin:$PATH\n
      setenv PYTHONPATH /workspace/hshu/anaconda3/bin:$PATH\n
      setenv MKL_NUM_THREADS 1\n"
      
      script[4] = sprintf("python simulation2_OnPLS.py %d",j)
      
      write.table(script,file = sprintf("/scratch/hshu/paper1_revision/setup2/jobs/job_simulation2_OnPLS%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("qsub /scratch/hshu/paper1_revision/setup2/jobs/job_simulation2_OnPLS%d.pbs",j),intern=F)
    }    
    
    
