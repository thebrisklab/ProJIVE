label=load('PAM50_label.txt');% subtype labels
X1_noisy = load('Expression.txt');%mRNA
X2_noisy = load('Methylation_above90.txt');

[p1,n] = size(X1_noisy);
p2 = size(X2_noisy,1);

w = floor(0.5*sqrt(n));

rng(0);
numfactors_nonstd_log(X1_noisy',10,floor(p1/2),floor(p1/40),1,'p2',1000,w,w);
rng(0);
numfactors_nonstd_log(X2_noisy',10,floor(p2/2),floor(p2/40),1,'p2',1000,w,w);
rng(0) 
numfactors_nonstd_log([X1_noisy' X2_noisy'],10,floor((p1+p2)/2),floor((p1+p2)/20),1,'p2',1000,w,w);

q1=3;
q2=3;
q12=4;

rng(0);
[chi_xy1, chi_y1, phi_y1, psi_y1, nu_y1, xi_xy1,   chi_xy2, chi_y2, phi_y2, psi_y2, nu_y2, xi_xy2] = kestimate_2_blocks(X1_noisy',X2_noisy',q12,q1,q2,0);

save(['chi_xy1_above90.txt'],'chi_xy1','-ascii');
save(['chi_y1_above90.txt'],'chi_y1','-ascii');
save(['phi_y1_above90.txt'],'phi_y1','-ascii');
save(['psi_y1_above90.txt'],'psi_y1','-ascii');
save(['nu_y1_above90.txt'],'nu_y1','-ascii');
save(['xi_xy1_above90.txt'],'xi_xy1','-ascii');

save(['chi_xy2_above90.txt'],'chi_xy2','-ascii');
save(['chi_y2_above90.txt'],'chi_y2','-ascii');
save(['phi_y2_above90.txt'],'phi_y2','-ascii');
save(['psi_y2_above90.txt'],'psi_y2','-ascii');
save(['nu_y2_above90.txt'],'nu_y2','-ascii');
save(['xi_xy2_above90.txt'],'xi_xy2','-ascii');




%%Block gdfm
xi_y1= nu_y1+xi_xy1;
xi_y2= nu_y2+xi_xy2;
psi_nu_y1= psi_y1+nu_y1;
psi_nu_y2= psi_y2+nu_y2;

%var_chi_xy1=round( norm(chi_xy1,'fro')^2/norm(chi_xy1,'fro')^2 ,3)
%var_chi_xy2=round( norm(chi_xy2,'fro')^2/norm(chi_xy2,'fro')^2 ,3)

var_chi_y1=round( norm(chi_y1,'fro')^2/norm(chi_xy1,'fro')^2 ,3)
var_chi_y2=round( norm(chi_y2,'fro')^2/norm(chi_xy2,'fro')^2 ,3)

var_phi_y1=round( norm(phi_y1,'fro')^2/norm(chi_xy1,'fro')^2 ,3)
var_phi_y2=round( norm(phi_y2,'fro')^2/norm(chi_xy2,'fro')^2 ,3)

var_psi_y1=round( norm(psi_y1,'fro')^2/norm(chi_xy1,'fro')^2 ,3)
var_psi_y2=round( norm(psi_y2,'fro')^2/norm(chi_xy2,'fro')^2 ,3)

var_nu_y1=round( norm(nu_y1,'fro')^2/norm(chi_xy1,'fro')^2 ,3)
var_nu_y2=round( norm(nu_y2,'fro')^2/norm(chi_xy2,'fro')^2 ,3)

var_xi_xy1=round( norm(xi_xy1,'fro')^2/norm(chi_xy1,'fro')^2 ,3)
var_xi_xy2=round( norm(xi_xy2,'fro')^2/norm(chi_xy2,'fro')^2 ,3)

var_xi_y1=round( norm(xi_y1,'fro')^2/norm(chi_xy1,'fro')^2 ,3)
var_xi_y2=round( norm(xi_y2,'fro')^2/norm(chi_xy2,'fro')^2 ,3)

var_psi_nu_y1=round( norm(psi_nu_y1,'fro')^2/norm(chi_xy1,'fro')^2 ,3)
var_psi_nu_y2=round( norm(psi_nu_y2,'fro')^2/norm(chi_xy2,'fro')^2 ,3)

rank(chi_xy1)%,1e-6)
rank(chi_xy2)%,1e-5)

rank(chi_y1)%,1e-6)
rank(chi_y2)%,1e-6)

rank(phi_y1)%,1e-6)
rank(phi_y2)%,1e-6)

rank(psi_y1)%,1e-6)
rank(psi_y2)%,1e-6)

rank(nu_y1)%,1e-6)
rank(nu_y2)%,1e-6)

rank(xi_xy1)%,1e-6)
rank(xi_xy2)%,1e-6)

rank(xi_y1)%,1e-6)
rank(xi_y2)%,1e-6)

rank(psi_nu_y1)%,1e-6)
rank(psi_nu_y2)%,1e-6)

round(SWISS(chi_xy1',[],label,1),3)
round(SWISS(chi_xy2',[],label,1),3)

round(SWISS(chi_y1',[],label,1),3)
round(SWISS(chi_y2',[],label,1),3)

round(SWISS(phi_y1',[],label,1),3)
round(SWISS(phi_y2',[],label,1),3)

round(SWISS(psi_y1',[],label,1),3)
round(SWISS(psi_y2',[],label,1),3)

round(SWISS(nu_y1',[],label,1),3)
round(SWISS(nu_y2',[],label,1),3)

round(SWISS(xi_xy1',[],label,1),3)
round(SWISS(xi_xy2',[],label,1),3)

round(SWISS(xi_y1',[],label,1),3)
round(SWISS(xi_y2',[],label,1),3)

round(SWISS(psi_nu_y1',[],label,1),3)
round(SWISS(psi_nu_y2',[],label,1),3)
  
