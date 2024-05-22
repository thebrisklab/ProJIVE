function simulation1_COBE(seed)
seed=num2str(seed);
Y_1=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle45_nstd1_Y_1.txt']);
Y_2=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle45_nstd1_Y_2.txt']);

[p_1,nn]=size(Y_1);
[p_2,nn]=size(Y_2);


%% COBE
rng(0);
addpath('/scratch/hshu/paper1_revision/demo_CIFA');

norm_Y_1=norm(Y_1,'fro');
Y_1_smat=Y_1/norm_Y_1;

norm_Y_2=norm(Y_2,'fro');
Y_2_smat=Y_2/norm_Y_2;

[U1 D1 V1]=svd(Y_1_smat,'econ');
[U2 D2 V2]=svd(Y_2_smat,'econ');

%denoise 
rank1=3;
rank2=3;

Y_1_smat_denoise=U1(:,1:rank1)*D1(1:rank1,1:rank1)*V1(:,1:rank1)';
Y_2_smat_denoise=U2(:,1:rank2)*D2(1:rank2,1:rank2)*V2(:,1:rank2)';

X_1_cobe=Y_1_smat_denoise*norm_Y_1;
X_2_cobe=Y_2_smat_denoise*norm_Y_2;

Y{1}=Y_1_smat_denoise';
Y{2}=Y_2_smat_denoise';

Ac=cobe(Y);%No common basis found


if sum(size(Ac))==0
  C_1_cobe=zeros(p_1,nn);
  C_2_cobe=zeros(p_2,nn);
else
  C_1_cobe=X_1_cobe*Ac*Ac';
  C_2_cobe=X_2_cobe*Ac*Ac';
end

D_1_cobe=X_1_cobe-C_1_cobe;
D_2_cobe=X_2_cobe-C_2_cobe;

X_1=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle45_nstd1_X_1.txt']);
X_2=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle45_nstd1_X_2.txt']);

C_1=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle45_nstd1_C_1.txt']);
C_2=load(['/scratch/hshu/paper1_revision/setup1/datafile/siml1_seed',seed,'_p900_angle45_nstd1_C_2.txt']);

D_1=X_1-C_1;
D_2=X_2-C_2;


errF_X_1=norm(X_1-X_1_cobe,'fro')/norm(X_1,'fro');
errF_X_2=norm(X_2-X_2_cobe,'fro')/norm(X_2,'fro');

errF_C_1=norm(C_1-C_1_cobe,'fro')/norm(C_1,'fro');
errF_C_2=norm(C_2-C_2_cobe,'fro')/norm(C_2,'fro');

errF_D_1=norm(D_1-D_1_cobe,'fro')/norm(D_1,'fro');
errF_D_2=norm(D_2-D_2_cobe,'fro')/norm(D_2,'fro');

err2_X_1=norm(X_1-X_1_cobe,2)/norm(X_1,2);
err2_X_2=norm(X_2-X_2_cobe,2)/norm(X_2,2);

err2_C_1=norm(C_1-C_1_cobe,2)/norm(C_1,2);
err2_C_2=norm(C_2-C_2_cobe,2)/norm(C_2,2);

err2_D_1=norm(D_1-D_1_cobe,2)/norm(D_1,2);
err2_D_2=norm(D_2-D_2_cobe,2)/norm(D_2,2);

output=[errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2,rank(C_1_cobe),rank(C_2_cobe)];

save(['siml1_seed' seed '_p900_angle45_nstd1_COBE.txt'],'output','-ascii');

end