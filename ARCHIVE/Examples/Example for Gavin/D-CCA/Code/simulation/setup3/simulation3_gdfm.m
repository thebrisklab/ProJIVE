  Y_1=load(['Y_1.txt']);
  Y_2=load(['Y_2.txt']);

  [p1,n] = size(Y_1);
  p2 = size(Y_2,1);
  
  w = floor(0.5*sqrt(n));
  
  q1=3;
  q2=5;
  q12=8;

  
  rng(0);
  [chi_xy1, chi_y1, phi_y1, psi_y1, nu_y1, xi_xy1,   chi_xy2, chi_y2, phi_y2, psi_y2, nu_y2, xi_xy2] = kestimate_2_blocks(Y_1',Y_2',q12,q1,q2,0);
  
  xi_y1= nu_y1+xi_xy1;
  xi_y2= nu_y2+xi_xy2;
  
  save(['chi_xy1_gdfm.txt'],'chi_xy1','-ascii');
  save(['chi_xy2_gdfm.txt'],'chi_xy2','-ascii');
  
  save(['chi_y1_gdfm.txt'],'chi_y1','-ascii');
  save(['chi_y2_gdfm.txt'],'chi_y2','-ascii');
  
  save(['phi_y1_gdfm.txt'],'phi_y1','-ascii');
  save(['phi_y2_gdfm.txt'],'phi_y2','-ascii');
  
  save(['psi_y1_gdfm.txt'],'psi_y1','-ascii');
  save(['psi_y2_gdfm.txt'],'psi_y2','-ascii');
  
  save(['nu_y1_gdfm.txt'],'nu_y1','-ascii');
  save(['nu_y2_gdfm.txt'],'nu_y2','-ascii');  
  
  save(['xi_y1_gdfm.txt'],'xi_y1','-ascii');
  save(['xi_y2_gdfm.txt'],'xi_y2','-ascii');
  
  save(['xi_xy1_gdfm.txt'],'xi_xy1','-ascii');
  save(['xi_xy2_gdfm.txt'],'xi_xy2','-ascii');

  
  
  
  %test correlations between psi_y1 and psi_y2
  t=1;
  for i=1:p1
    for j=1:p2
      nlogP(t)=-log(test_zero_corr(psi_y1(:,i),psi_y2(:,j)));
      t=t+1;
    end
  end
  
  save(['siml3_gdfm_nlogP_1.txt'],'nlogP','-ascii');
  
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
  
  save(['siml3_gdfm_nlogP_2.txt'],'nlogP','-ascii');


