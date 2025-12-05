%analyze BRCA data 
rng(0)
label=load('PAM50_label.txt');% subtype labels
num_1=sum(label==1)%112
num_2=sum(label==2)%331
num_3=sum(label==3)%162
num_4=sum(label==4)%55
%%%%

X1_noisy=load('Expression.txt');%mRNA
X2_noisy=load('Methylation_below90.txt');


size(X1_noisy)
size(X2_noisy)

[p1,nn]=size(X1_noisy)
[p2,~]=size(X2_noisy)

%% proposed method
X1_ours=load('X_1_hat_ours_below90.txt');
X2_ours=load('X_2_hat_ours_below90.txt');
C1_ours=load('C_1_hat_ours_below90.txt');
C2_ours=load('C_2_hat_ours_below90.txt');
D1_ours=load('D_1_hat_ours_below90.txt');
D2_ours=load('D_2_hat_ours_below90.txt');


%%% explained variation 
var_C1_ours=round( norm(C1_ours,'fro')^2/norm(X1_ours,'fro')^2  ,3)
var_C2_ours=round( norm(C2_ours,'fro')^2/norm(X2_ours,'fro')^2  ,3)

var_D1_ours=round( norm(D1_ours,'fro')^2/norm(X1_ours,'fro')^2  ,3)
var_D2_ours=round( norm(D2_ours,'fro')^2/norm(X2_ours,'fro')^2  ,3)

%%% rank of matrix
rank(X1_ours)
rank(X2_ours)
rank(C1_ours)
rank(C2_ours)
rank(D1_ours)
rank(D2_ours)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SWISS
round(SWISS(X1_noisy,[],label,1),3) 
round(SWISS(X2_noisy,[],label,1),3)

round(SWISS(X1_ours,[],label,1),3) 
round(SWISS(X2_ours,[],label,1),3) 

round(SWISS(C1_ours,[],label,1),3) 
round(SWISS(C2_ours,[],label,1),3) 

round(SWISS(D1_ours,[],label,1),3) 
round(SWISS(D2_ours,[],label,1),3) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   other methods
%% JIVE without the orthongonal constraint on distinct matrices (jive)
%See jive_analyze_below90.R
X1_jive=load('X_1_jive_below90.txt');
X2_jive=load('X_2_jive_below90.txt');

C1_jive=load('C_1_jive_below90.txt');
C2_jive=load('C_2_jive_below90.txt');

D1_jive=load('D_1_jive_below90.txt');
D2_jive=load('D_2_jive_below90.txt');

%matrix rank
rank(X1_jive)
rank(X2_jive)
rank(C1_jive)
rank(C2_jive)
rank(D1_jive)
rank(D2_jive)

%SWISS
round( SWISS(X1_jive,[],label,1) ,3)
round( SWISS(X2_jive,[],label,1) ,3)

round( SWISS(C1_jive,[],label,1) ,3)
round( SWISS(C2_jive,[],label,1) ,3)

round( SWISS(D1_jive,[],label,1) ,3)
round( SWISS(D2_jive,[],label,1) ,3)

% explained variation 
var_C1_jive=round( norm(C1_jive,'fro')^2/norm(X1_jive,'fro')^2   ,3)
var_C2_jive=round( norm(C2_jive,'fro')^2/norm(X2_jive,'fro')^2   ,3)

var_D1_jive=round( norm(D1_jive,'fro')^2/norm(X1_jive,'fro')^2   ,3)
var_D2_jive=round( norm(D2_jive,'fro')^2/norm(X2_jive,'fro')^2   ,3)

%% JIVE with the orthongonal constraint on distinct matrices (rJIVE)
%See jive_analyze_below90.R
X1_rJIVE=load('X_1_rjive_below90.txt');
X2_rJIVE=load('X_2_rjive_below90.txt');

C1_rJIVE=load('C_1_rjive_below90.txt');
C2_rJIVE=load('C_2_rjive_below90.txt');

D1_rJIVE=load('D_1_rjive_below90.txt');
D2_rJIVE=load('D_2_rjive_below90.txt');

%matrix rank
rank(X1_rJIVE)
rank(X2_rJIVE)
rank(C1_rJIVE)
rank(C2_rJIVE)
rank(D1_rJIVE)
rank(D2_rJIVE)

%SWISS
round( SWISS(X1_rJIVE,[],label,1) ,3)
round( SWISS(X2_rJIVE,[],label,1) ,3)

round( SWISS(C1_rJIVE,[],label,1) ,3)
round( SWISS(C2_rJIVE,[],label,1) ,3)

round( SWISS(D1_rJIVE,[],label,1) ,3)
round( SWISS(D2_rJIVE,[],label,1) ,3)

% explained variation 
var_C1_rJIVE=round( norm(C1_rJIVE,'fro')^2/norm(X1_rJIVE,'fro')^2   ,3)
var_C2_rJIVE=round( norm(C2_rJIVE,'fro')^2/norm(X2_rJIVE,'fro')^2   ,3)

var_D1_rJIVE=round( norm(D1_rJIVE,'fro')^2/norm(X1_rJIVE,'fro')^2   ,3)
var_D2_rJIVE=round( norm(D2_rJIVE,'fro')^2/norm(X2_rJIVE,'fro')^2   ,3)




%% ajive
rng(0)
addpath('/rsrch2/biostatistics/hshu/paper1/AJIVECode');

vecr = [2, 3];

datablock{1} = X1_noisy;
datablock{2} = X2_noisy; 

paramstruct = struct('imean',[0 0],'ioutput', [0,0,0,0,0,0, 1, 1, 0],'iplot',[0,0]);

outstruct = AJIVEMainMJ(datablock, vecr,paramstruct);


C1_ajive = outstruct.MatrixJoint{1};
C2_ajive = outstruct.MatrixJoint{2};
D1_ajive = outstruct.MatrixIndiv{1};
D2_ajive= outstruct.MatrixIndiv{2};

X1_ajive=C1_ajive+D1_ajive;
X2_ajive=C2_ajive+D2_ajive;

save(['C_1_hat_AJIVE_below90.txt'],'C1_ajive','-ascii');
save(['C_2_hat_AJIVE_below90.txt'],'C2_ajive','-ascii');

save(['D_1_hat_AJIVE_below90.txt'],'D1_ajive','-ascii');
save(['D_2_hat_AJIVE_below90.txt'],'D2_ajive','-ascii');

save(['X_1_hat_AJIVE_below90.txt'],'X1_ajive','-ascii');
save(['X_2_hat_AJIVE_below90.txt'],'X2_ajive','-ascii');

%matrix rank
rank(X1_ajive)
rank(X2_ajive)

rank(C1_ajive)
rank(C2_ajive)

rank(D1_ajive)
rank(D2_ajive)
% explained variation 
var_C1_ajive=round( norm(C1_ajive,'fro')^2/norm(X1_ajive,'fro')^2  ,3)
var_C2_ajive=round( norm(C2_ajive,'fro')^2/norm(X2_ajive,'fro')^2  ,3)

var_D1_ajive=round( norm(D1_ajive,'fro')^2/norm(X1_ajive,'fro')^2  ,3)
var_D2_ajive=round( norm(D2_ajive,'fro')^2/norm(X2_ajive,'fro')^2  ,3)


%SWISS score

round( SWISS(X1_ajive,[],label,1) ,3)
round( SWISS(X2_ajive,[],label,1) ,3)

round( SWISS(C1_ajive,[],label,1) ,3)
round( SWISS(C2_ajive,[],label,1) ,3)

round( SWISS(D1_ajive,[],label,1) ,3)
round( SWISS(D2_ajive,[],label,1) ,3)


%% COBE
addpath('/rsrch2/biostatistics/hshu/software/demo_CIFA');
%To obtain the clean data, we use SORTE+standard PCA as suggested in Remark 4 of the COBE paper
rng(0)

r1_cobe=SORTE(X1_noisy)%2
r2_cobe=SORTE(X2_noisy)%1

[U1 D1 V1]=svd(X1_noisy,'econ');
[U2 D2 V2]=svd(X2_noisy,'econ');

X1_cobe=U1(:,1:r1_cobe)*D1(1:r1_cobe,1:r1_cobe)*V1(:,1:r1_cobe)';
X2_cobe=U2(:,1:r2_cobe)*D2(1:r2_cobe,1:r2_cobe)*V2(:,1:r2_cobe)';

rank(X1_cobe)
rank(X2_cobe)

X1_cobe_norm=norm(X1_cobe,'fro');
X2_cobe_norm=norm(X2_cobe,'fro');

X1_cobe_s=X1_cobe/X1_cobe_norm;
X2_cobe_s=X2_cobe/X2_cobe_norm;

Y{1}=X1_cobe_s';
Y{2}=X2_cobe_s';

Ac=cobe(Y);%No common basis found.

round(SWISS(X1_cobe,[],label,1),3) %
round(SWISS(X2_cobe,[],label,1),3) %

%% OnPLS
X1_OnPLS=load('X_1_hat_OnPLS_below90.txt');
X2_OnPLS=load('X_2_hat_OnPLS_below90.txt');

C1_OnPLS=load('C_1_hat_OnPLS_below90.txt');
C2_OnPLS=load('C_2_hat_OnPLS_below90.txt');

D1_OnPLS=load('D_1_hat_OnPLS_below90.txt');
D2_OnPLS=load('D_2_hat_OnPLS_below90.txt');

%matrix rank
rank(X1_OnPLS)
rank(X2_OnPLS)
rank(C1_OnPLS)
rank(C2_OnPLS)
rank(D1_OnPLS)
rank(D2_OnPLS)

% explained variation 
var_C1_OnPLS=round( norm(C1_OnPLS,'fro')^2/norm(X1_OnPLS,'fro')^2 ,3)
var_C2_OnPLS=round( norm(C2_OnPLS,'fro')^2/norm(X2_OnPLS,'fro')^2 ,3)

var_D1_OnPLS=round( norm(D1_OnPLS,'fro')^2/norm(X1_OnPLS,'fro')^2 ,3)
var_D2_OnPLS=round( norm(D2_OnPLS,'fro')^2/norm(X2_OnPLS,'fro')^2 ,3)
%SWISS score

round(SWISS(X1_OnPLS,[],label,1),3) %
round(SWISS(X2_OnPLS,[],label,1),3) %

round(SWISS(C1_OnPLS,[],label,1),3) %
round(SWISS(C2_OnPLS,[],label,1),3) %

round(SWISS(D1_OnPLS,[],label,1),3) %
round(SWISS(D2_OnPLS,[],label,1),3) %


%% DISCO-SCA
norm_X1_noisy=norm(X1_noisy,'fro');
norm_X2_noisy=norm(X2_noisy,'fro');

  
X1_noisy_bscaled=X1_noisy/norm_X1_noisy;
X2_noisy_bscaled=X2_noisy/norm_X2_noisy;

  
con_data=[X1_noisy_bscaled',X2_noisy_bscaled'];
row1=[ones(1,p1),2*ones(1,p2)];
con_data=[row1;con_data];
size(con_data)
save(['concatenated_data_blockscaled_below90.txt'],'con_data','-ascii');

addpath('E:\Hai Windows\work\softwares\matlab_toolbox\DISCO-SCA_Source');
%%%First to decide the number of total common and distinct components
rng(0)
DISCO_SCA
%%%We choose total number=13 by DISCOrankselectionbelow90_ScreePlots.png
%%%Next, we run the DISCO-SCA
rng(0)
DISCO_SCA
%see the parameter settings in DISCOsettings_below90.JPG
%we choose VARIMAX instead of EQUAMAX because the latter one does not
%converge when max # of iteration is set to be 50000.

%choose 7 distinctive components according to
%DISCOresultbelow90_status_components.png
result=load('DISCOresultbelow90_7.mat')
result.DISCO_SCA.AnalysisOptions.Target
%     1     1     1     1     0     0     0     1     1     1     1     1     1
%     0     0     0     0     1     1     1     1     1     1     1     1     1
%The 8~13-th columns are all 1, which indicates the 8~13-th components are common.

result.DISCO_SCA.AnalysisOutput

T=result.DISCO_SCA.AnalysisOutput.RotatedScores;
P=result.DISCO_SCA.AnalysisOutput.RotatedLoadings;


X1_DISCO=(T*P(1:p1,:)')'*norm_X1_noisy;
X2_DISCO=(T*P((p1+1):(p1+p2),:)')'*norm_X2_noisy;

C1_DISCO=(T(:,8:13)*P(1:p1,8:13)')'*norm_X1_noisy;
C2_DISCO=(T(:,8:13)*P((p1+1):(p1+p2),8:13)')'*norm_X2_noisy;

D1_DISCO=X1_DISCO-C1_DISCO;
D2_DISCO=X2_DISCO-C2_DISCO;


save(['C_1_hat_DISCO_below90.txt'],'C1_DISCO','-ascii');
save(['C_2_hat_DISCO_below90.txt'],'C2_DISCO','-ascii');

save(['D_1_hat_DISCO_below90.txt'],'D1_DISCO','-ascii');
save(['D_2_hat_DISCO_below90.txt'],'D2_DISCO','-ascii');

save(['X_1_hat_DISCO_below90.txt'],'X1_DISCO','-ascii');
save(['X_2_hat_DISCO_below90.txt'],'X2_DISCO','-ascii');

%matrix rank
rank(X1_DISCO)%13
rank(X2_DISCO)%13
rank(C1_DISCO)%6
rank(C2_DISCO)%6
rank(D1_DISCO)%7
rank(D2_DISCO)%7

% explained variation 
var_C1_DISCO=round( norm(C1_DISCO,'fro')^2/norm(X1_DISCO,'fro')^2 ,3)%0.732
var_C2_DISCO=round( norm(C2_DISCO,'fro')^2/norm(X2_DISCO,'fro')^2 ,3)%0.571

var_D1_DISCO=round( norm(D1_DISCO,'fro')^2/norm(X1_DISCO,'fro')^2 ,3)%0.268
var_D2_DISCO=round( norm(D2_DISCO,'fro')^2/norm(X2_DISCO,'fro')^2 ,3)%0.429
%SWISS score

round(SWISS(X1_DISCO,[],label,1),3) %0.526
round(SWISS(X2_DISCO,[],label,1),3) %0.663


round(SWISS(C1_DISCO,[],label,1),3) %0.447
round(SWISS(C2_DISCO,[],label,1),3) %0.400


round(SWISS(D1_DISCO,[],label,1),3) %0.935
round(SWISS(D2_DISCO,[],label,1),3) %0.992

