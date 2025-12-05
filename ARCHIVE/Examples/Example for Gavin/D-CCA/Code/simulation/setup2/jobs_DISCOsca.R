#nautilus
flist=  list.files("/rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2", pattern = "*_p900_angle45_nstd1_DISCO.txt")
fdone=as.integer(sub('.*siml2_seed*(.*?) *_p900_angle45_nstd1_DISCO.txt.*','\\1',flist))


for(j in 1:1000)
{
  
  script=NULL
  script[1]=sprintf("#PBS -N DISCO2")
  script[2] ="#PBS -l procs=1,mem=4gb\n
  #PBS -l walltime=120:00:00\n
  #PBS -d /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2\n
  #PBS -o /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2/output\n
  #PBS -e /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2/output\n"
  
  script[3] = "cat\n"  
  
  script[4] = "module load R/3.3.3-shlib\n" 
  
  script[5] = sprintf("Rscript simulation2_DISCOsca.R %d",j)
  
  write.table(script,file = sprintf("/rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2/jobs/job_simulation_disco%d.pbs",j),row.names=F,col.names = F,quote=F)
  
  system(sprintf("msub < /rsrch2/biostatistics/hshu/paper1/simulation1/setup2_paper_revision2/jobs/job_simulation_disco%d.pbs",j),intern=F)
  
  
  
}
