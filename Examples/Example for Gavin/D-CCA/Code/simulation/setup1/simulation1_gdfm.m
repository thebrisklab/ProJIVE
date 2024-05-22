function simulation1_gdfm(seed_i)

for seed_ii=(40*(seed_i-1)+1):(40*seed_i)
seed=num2str(seed_ii);

theta_list=[0,45];

for theta_i=1:2

  theta=num2str(theta_list(theta_i));
  
  Y_1=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle',theta,'_nstd1_Y_1.txt']);
  Y_2=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle',theta,'_nstd1_Y_2.txt']);
  
  X_1=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle',theta,'_nstd1_X_1.txt']);
  X_2=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle',theta,'_nstd1_X_2.txt']);
  
  
  [p1,n] = size(Y_1);
  p2 = size(Y_2,1);
  
  w = floor(0.5*sqrt(n));
  
  q1=3;
  q2=3;
  if theta_i==1
    q12=5;
  else
    q12=6;
  end
  
  rng(0);
  [chi_xy1, chi_y1, phi_y1, psi_y1, nu_y1, xi_xy1,   chi_xy2, chi_y2, phi_y2, psi_y2, nu_y2, xi_xy2] = kestimate_2_blocks(Y_1',Y_2',q12,q1,q2,0);
  
  xi_y1= nu_y1+xi_xy1;
  xi_y2= nu_y2+xi_xy2;
  
  
  errF_chi_xy1=norm(X_1-chi_xy1','fro')/norm(X_1,'fro');
  errF_chi_xy2=norm(X_2-chi_xy2','fro')/norm(X_2,'fro');
  
  err2_chi_xy1=norm(X_1-chi_xy1',2)/norm(X_1,2);
  err2_chi_xy2=norm(X_2-chi_xy2',2)/norm(X_2,2);
  
  
  errF_chi_y1=norm(X_1-chi_y1','fro')/norm(X_1,'fro');
  errF_chi_y2=norm(X_2-chi_y2','fro')/norm(X_2,'fro');
  
  err2_chi_y1=norm(X_1-chi_y1',2)/norm(X_1,2);
  err2_chi_y2=norm(X_2-chi_y2',2)/norm(X_2,2);
  
  
  normF_nu_y1_X=norm(nu_y1,'fro')/norm(X_1,'fro');
  normF_nu_y2_X=norm(nu_y2,'fro')/norm(X_2,'fro');
  
  norm2_nu_y1_X=norm(nu_y1,2)/norm(X_1,2);
  norm2_nu_y2_X=norm(nu_y2,2)/norm(X_2,2);
  
  
  normF_nu_y1_chi=norm(nu_y1,'fro')/norm(chi_xy1,'fro');
  normF_nu_y2_chi=norm(nu_y2,'fro')/norm(chi_xy2,'fro');
  
  norm2_nu_y1_chi=norm(nu_y1,2)/norm(chi_xy1,2);
  norm2_nu_y2_chi=norm(nu_y2,2)/norm(chi_xy2,2);
  
  
  output=[max(max(abs(phi_y1))), max(max(abs(phi_y2))), ...
                      errF_chi_xy1, errF_chi_xy2, err2_chi_xy1, err2_chi_xy2, ...
                      errF_chi_y1, errF_chi_y2, err2_chi_y1, err2_chi_y2, ...
                      normF_nu_y1_X, normF_nu_y2_X, norm2_nu_y1_X, norm2_nu_y2_X, ...
                      normF_nu_y1_chi, normF_nu_y2_chi, norm2_nu_y1_chi, norm2_nu_y2_chi];
  
  save(['siml1_seed', seed, '_p900_angle',theta,'_nstd1_gdfm.txt'],'output','-ascii');
  
  %test correlations between psi_y1 and psi_y2
  t=1;
  for i=1:p1
    for j=1:p2
      nlogP(t)=-log(test_zero_corr(psi_y1(:,i),psi_y2(:,j)));
      t=t+1;
    end
  end
  
  save(['siml1_seed', seed, '_p900_angle',theta,'_nstd1_gdfm_nlogP_1.txt'],'nlogP','-ascii');
  
  %test correlations between psi_y1+nu_y1 and psi_y2+nu_y2
  
  psi_nu_y1=psi_y1+nu_y1;
  psi_nu_y2=psi_y2+nu_y2;
  
  t=1;
  for i=1:p1
    for j=1:p2
      nlogP(t)=-log(test_zero_corr(psi_nu_y1(:,i),psi_nu_y2(:,j)));
      t=t+1;
    end
  end
  
  save(['siml1_seed', seed, '_p900_angle',theta,'_nstd1_gdfm_nlogP_2.txt'],'nlogP','-ascii');

end

end