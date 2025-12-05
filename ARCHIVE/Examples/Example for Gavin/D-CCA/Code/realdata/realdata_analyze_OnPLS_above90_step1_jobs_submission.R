for(i in 1:10){
  for(j in 0:10){
    for(k in 0:10){
      script=NULL
      script[1]=sprintf("#PBS -N OnPLSa")
      script[2] ="#PBS -l procs=1,pmem=2gb\n
      #PBS -l walltime=2:00:00\n"    

      script[3] = "cd /scratch/hshu/paper1_revision/realdata_revision2\n
      setenv PATH /workspace/hshu/anaconda3/bin:$PATH\n
      setenv PYTHONPATH /workspace/hshu/anaconda3/bin:$PATH\n
      setenv MKL_NUM_THREADS 1\n"
      
      script[4] = sprintf("python realdata_analyze_OnPLS_above90_step1.py %d %d %d",i,j,k)
      
      write.table(script,file = sprintf("/scratch/hshu/paper1_revision/realdata_revision2/jobs/job_above90_OnPLS%d_%d_%d.pbs",i,j,k),row.names=F,col.names = F,quote=F)
      
      system(sprintf("qsub /scratch/hshu/paper1_revision/realdata_revision2/jobs/job_above90_OnPLS%d_%d_%d.pbs",i,j,k),intern=F)
      
    }
  }
}

  


