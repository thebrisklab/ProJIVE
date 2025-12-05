function [pvalue]=test_zero_corr(x,y,tail)
%assume x,y are column vectors
if nargin == 2
  tail='two';
end

x=x-mean(x);
y=y-mean(y);

if min(max(abs(x)),max(abs(y)))<1e-8
    pvalue=1;

else
  
  if any(strcmp(tail,{'two','left','right'}))
  
      n=length(x);
      x2=x.^2;
      y2=y.^2;
      sum_x2=sum(x2);
      sum_y2=sum(y2);
      rho=x'*y/sqrt(sum_x2*sum_y2);
      tau=sqrt(    n* x2'*y2/(sum_x2*sum_y2)  );
      test_stat=sqrt(n)*rho/tau;
      
      
      if strcmp(tail,'two')
        pvalue=2*normcdf(-abs(test_stat));
      elseif strcmp(tail,'left')
        pvalue=normcdf(test_stat);
      else
        pvalue=1-normcdf(test_stat);
      end


  else
    error('The argument "tail" must be one of ''two'', ''left'' and ''right'' .');
  end
  
end