  #dqshtc
flist=  list.files("/scratch/hshu/paper1_revision/setup2/revision2", pattern = "*_p900_angle45_nstd1_gdfm_nlogP_2.txt")
fdone=as.integer(sub('.*siml2_seed*(.*?) *_p900_angle45_nstd1_gdfm_nlogP_2.txt.*','\\1',flist))#933

#    for(j in 1:25)
    for(j in setdiff(1:1000,fdone))
    {
      script=NULL
      script[1]=sprintf("#PBS -N gdfm2")
      script[2] ="#PBS -l procs=1,pmem=2gb\n
#PBS -l walltime=24:00:00"      
      
      script[3] = "cd /scratch/hshu/paper1_revision/setup2/revision2\n
      module load matlab/2014b\n"
      
      script[4] = sprintf("matlab -nodisplay -nodesktop -singleCompThread -r 'simulation2_gdfm(%d)'",j)
      
      write.table(script,file = sprintf("/scratch/hshu/paper1_revision/setup2/jobs/job_simulation2_gdfm%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("qsub /scratch/hshu/paper1_revision/setup2/jobs/job_simulation2_gdfm%d.pbs",j),intern=F)
    }    
    
    
    
    
    
    
    