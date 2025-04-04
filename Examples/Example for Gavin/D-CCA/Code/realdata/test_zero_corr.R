test_zero_corr<-function(x,y,tail='two'){
  x=x-mean(x)
  y=y-mean(y)
  if( min(max(abs(x)),max(abs(y)))<1e-8 ){
    pvalue=1
    return(pvalue)
  }else{
    if(tail%in%c('two','left','right'))
    {
      n=length(x)
      x.2=x^2
      y.2=y^2
      sum.x.2=sum(x.2)
      sum.y.2=sum(y.2)
      rho=t(x)%*%y/sqrt(sum.x.2*sum.y.2)
      tau=sqrt(    n* t(x.2)%*%y.2/(sum.x.2*sum.y.2)  )
      test_stat=sqrt(n)*rho/tau
      if(tail=='two'){
        pvalue=2*pnorm(-abs(test_stat))
      }else if(tail=='left'){
        pvalue=pnorm(test_stat)
      }else{
        pvalue=1-pnorm(test_stat)
      }
      return(pvalue)
    }else{
      stop("'The argument \"tail\" must be one of \"two\", \"left\" and \"right\" .'")
    }
  }
  
}



