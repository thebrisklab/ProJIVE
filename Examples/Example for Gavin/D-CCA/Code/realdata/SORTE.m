function rank_signal = SORTE(X, k_start, prop)
    % See Section 4.3 in He, Z., Cichocki, A., Xie, S. and Choi, K., 2010. Detecting the number of clusters in n-way probabilistic clustering. 
    % IEEE Transactions on Pattern Analysis and Machine Intelligence, 32(11), pp.2006-2021.
    % prop: see the "p" in Remark 3 of the paper
    
    if nargin < 2
        k_start=1;
    end
    if nargin < 3
        prop=0.99;
    end
    
    T=size(X,2);
    lambda = (svd(X).^2)/T;
    
    n_lambda=length(lambda);
    if n_lambda<T
        lambda=[lambda;zeros(T-n_lambda,1)];
    else
        lambda=lambda(1:T);
    end
    
    
    sum_lambda=sum(lambda);
    sum_lambda_prop=sum_lambda*prop;
    
    J=T-2;
    psum_lambda=sum_lambda;
    while psum_lambda>sum_lambda_prop
       J=J-1;
       psum_lambda=sum_lambda-sum(lambda(J:T));
    end
    
    
    
    delta=lambda(1:(T-1))-lambda(2:T);
    sigma2=zeros(1,J+1);
    SORTE_k=zeros(1,J);
    
    sigma2(1)=var(delta,1);
    for k=1:J
        sigma2(k+1)=var(delta((k+1):end),1);
        SORTE_k(k)=sigma2(k+1)/sigma2(k);
    end

    
    [~,rank_signal]=min(SORTE_k(k_start:end));
    rank_signal=rank_signal+k_start-1;
end
