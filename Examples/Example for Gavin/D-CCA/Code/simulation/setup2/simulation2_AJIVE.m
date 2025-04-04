function simulation2_AJIVE(seed)
seed=num2str(seed);

Y_1=load(['/scratch/hshu/paper1_revision/setup2/datafile/siml2_seed',seed,'_p900_angle45_nstd1_Y_1.txt']);
Y_2=load(['/scratch/hshu/paper1_revision/setup2/datafile/siml2_seed',seed,'_p900_angle45_nstd1_Y_2.txt']);

rng(0);
addpath('/scratch/hshu/paper1_revision/AJIVECode');

vecr = [3, 5];

datablock{1} = Y_1;
datablock{2} = Y_2; 

paramstruct = struct('imean',[0 0],'ioutput', [0,0,0,0,0,0, 1, 1, 0],'iplot',[0,0]);

outstruct = AJIVEMainMJ(datablock, vecr,paramstruct);


C_1_ajive = outstruct.MatrixJoint{1};
C_2_ajive = outstruct.MatrixJoint{2};
D_1_ajive = outstruct.MatrixIndiv{1};
D_2_ajive= outstruct.MatrixIndiv{2};

X_1_ajive=C_1_ajive+D_1_ajive;
X_2_ajive=C_2_ajive+D_2_ajive;

X_1=load(['/scratch/hshu/paper1_revision/setup2/datafile/siml2_seed',seed,'_p900_angle45_nstd1_X_1.txt']);
X_2=load(['/scratch/hshu/paper1_revision/setup2/datafile/siml2_seed',seed,'_p900_angle45_nstd1_X_2.txt']);

C_1=load(['/scratch/hshu/paper1_revision/setup2/datafile/siml2_seed',seed,'_p900_angle45_nstd1_C_1.txt']);
C_2=load(['/scratch/hshu/paper1_revision/setup2/datafile/siml2_seed',seed,'_p900_angle45_nstd1_C_2.txt']);

D_1=X_1-C_1;
D_2=X_2-C_2;


errF_X_1=norm(X_1-X_1_ajive,'fro')/norm(X_1,'fro');
errF_X_2=norm(X_2-X_2_ajive,'fro')/norm(X_2,'fro');

errF_C_1=norm(C_1-C_1_ajive,'fro')/norm(C_1,'fro');
errF_C_2=norm(C_2-C_2_ajive,'fro')/norm(C_2,'fro');

errF_D_1=norm(D_1-D_1_ajive,'fro')/norm(D_1,'fro');
errF_D_2=norm(D_2-D_2_ajive,'fro')/norm(D_2,'fro');

err2_X_1=norm(X_1-X_1_ajive,2)/norm(X_1,2);
err2_X_2=norm(X_2-X_2_ajive,2)/norm(X_2,2);

err2_C_1=norm(C_1-C_1_ajive,2)/norm(C_1,2);
err2_C_2=norm(C_2-C_2_ajive,2)/norm(C_2,2);

err2_D_1=norm(D_1-D_1_ajive,2)/norm(D_1,2);
err2_D_2=norm(D_2-D_2_ajive,2)/norm(D_2,2);

output=[errF_X_1, err2_X_1, errF_X_2, err2_X_2, errF_C_1, err2_C_1, errF_C_2, err2_C_2, errF_D_1, err2_D_1, errF_D_2, err2_D_2,rank(C_1_ajive),rank(C_2_ajive)];

save(['siml2_seed' seed '_p900_angle45_nstd1_AJIVE.txt'],'output','-ascii');

end
