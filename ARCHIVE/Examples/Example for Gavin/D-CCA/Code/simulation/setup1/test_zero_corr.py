from numpy.linalg import norm
import numpy as np
import scipy.stats as stats

def test_zero_corr(x, y, tail='two'):  
    #The studentized test for zero correlation based on the normal approximation 
    #which is proposed by DiCiccio, C. J. and Romano, J. P. (2017), “Robust Permutation Tests For Correlation And Regression Coefficients,” JASA.
    
    '''
    Assume x, y are mean zero.
    '''
    
    x = x - np.mean(x)
    y = y - np.mean(y)
    
    if np.min([np.max(abs(x)),np.max(abs(y))])<1e-8:
        pvalue = 1 # one of x, y is a constant
        return pvalue
    else:
        if tail in ['two', 'left', 'right']:
            n = len(x)
            rho = x.T @ y / ( norm(x) * norm(y) )
            tau = ( n*(x**2).T @ y**2 / (sum(x**2) * sum(y**2)) )**0.5
            test_stat = n**0.5 * rho / tau
            
            if tail=='two':
                pvalue = 2 * stats.norm.cdf(-abs(test_stat))                
            elif tail=='left':
                pvalue = stats.norm.cdf(test_stat)
            else:
                pvalue = 1 - stats.norm.cdf(test_stat)
                           
            return pvalue
            
        else:
            raise ValueError('The argument "tail" must be one of "two", "left" and "right" .') 
      