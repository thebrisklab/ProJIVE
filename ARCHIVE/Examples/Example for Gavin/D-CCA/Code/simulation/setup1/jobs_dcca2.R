#python3/3.5.1-0-anaconda
    for(j in 1:1000)
    {
      script=NULL
      script[1]=sprintf("#BSUB -J dcca1_2")
      script[2] ="#BSUB -W 23:59\n
#BSUB -q medium\n
#BSUB -n 1\n
#BSUB -M 8192\n
#BSUB -R rusage[mem=8192]\n
#BSUB -N \n"      
      
      script[3] = "#BSUB -o /rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/output\n
#BSUB -e /rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/output\n
#BSUB -cwd /rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1\n"
      
      script[4] = sprintf("python simulation1_dcca2.py %d",j)
      
      write.table(script,file = sprintf("/rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/jobs/job_simulation1%d.pbs",j),row.names=F,col.names = F,quote=F)
      
      system(sprintf("bsub < /rsrch2/biostatistics/hshu/paper1/simulation1/setup1_paper_revision1/jobs/job_simulation1%d.pbs",j),intern=F)
    }
