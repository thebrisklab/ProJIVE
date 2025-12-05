Y_1=load('Y_1.txt');
Y_2=load('Y_2.txt');
X_1=load('X_1.txt');
X_2=load('X_2.txt');
C_1=load('C_1.txt');
C_2=load('C_2.txt');
D_1=load('D_1.txt');
D_2=load('D_2.txt');
E_1=load('E_1.txt');
E_2=load('E_2.txt');



scale=max(max(max(abs(X_1))),max(max(abs(X_2)))  );


[ha, pos] =tight_subplot(2,6,[.03 .03],[.05 .02],[.05 .02]) ;
%subplot(2,6,1)
axes(ha(1)); 
      imagesc(Y_1)
 colormap(darkb2r(  -scale,scale   ))
  set(gca,'XTick',[]);
%colorbar

%subplot(2,6,2)
axes(ha(2)); 
      imagesc(X_1)
 colormap(darkb2r(  -scale,scale   ))
  set(gca,'XTick',[],'YTick',[]);
%colorbar
 
%subplot(2,6,3)
axes(ha(3)); 
   imagesc(C_1)
 colormap(darkb2r(-scale,scale ))
  set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,6,4)
axes(ha(4)); 
      imagesc(D_1)
 colormap(darkb2r(  -scale,scale   ))
  set(gca,'XTick',[],'YTick',[]);
%colorbar
    
%subplot(2,6,5)
axes(ha(5)); 
      imagesc(E_1)
 colormap(darkb2r(  -scale,scale   ))
  set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,6,11)
 axes(ha(6)); 
 %leave it blank 
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,6,6)
axes(ha(7)); 
      imagesc(Y_2)
 colormap(darkb2r(  -scale,scale   ))
 
%colorbar
    

%subplot(2,6,7)
axes(ha(8)); 
      imagesc(X_2)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,6,8)
axes(ha(9)); 
      imagesc(C_2)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar
    


%subplot(2,6,9)
axes(ha(10)); 
   imagesc(D_2)
 colormap(darkb2r(-scale,scale ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,6,10)
axes(ha(11)); 
      imagesc(E_2)
 colormap(darkb2r(  -scale,scale   ))
set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,6,12)
 axes(ha(12)); 
 %leave it blank   
 set(gca,'XTick',[],'YTick',[]);
%colorbar



%% proposed method
X_1_hat=load('X_1_hat_est_rank.txt');
X_2_hat=load('X_2_hat_est_rank.txt');
C_1_hat=load('C_1_hat_est_rank.txt');
C_2_hat=load('C_2_hat_est_rank.txt');
D_1_hat=load('D_1_hat_est_rank.txt');
D_2_hat=load('D_2_hat_est_rank.txt');


[ha, pos] =tight_subplot(2,6,[.03 .03],[.05 .02],[.05 .02]) ;
%subplot(2,5,1)
 axes(ha(1)); 
      imagesc(X_1_hat)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[]);
%colorbar
 
%subplot(2,5,2)
 axes(ha(2)); 
      imagesc(C_1_hat)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,5,3)
 axes(ha(3)); 
   imagesc(D_1_hat)
 colormap(darkb2r(-scale,scale ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,5,4)
 axes(ha(4)); 
%leave it blank
 set(gca,'XTick',[],'YTick',[]);
%colorbar
    
%subplot(2,5,5)
 axes(ha(5)); 
%leave it blank
 set(gca,'XTick',[],'YTick',[]);
%colorbar

%subplot(2,5,5)
 axes(ha(6)); 
%leave it blank
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,5,6)
 axes(ha(7)); 
      imagesc(X_2_hat)
 colormap(darkb2r(  -scale,scale   ))

%colorbar
    

%subplot(2,5,7)
 axes(ha(8)); 
      imagesc(C_2_hat)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar

%subplot(2,5,8)
 axes(ha(9)); 
      imagesc(D_2_hat)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar
    


%subplot(2,5,9)
 axes(ha(10)); 
 %leave it blank 
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,5,10)
 axes(ha(11)); 
 %leave it blank   
 set(gca,'XTick',[],'YTick',[]);
%colorbar

%subplot(2,5,10)
 axes(ha(12)); 
 %leave it blank   
 set(gca,'XTick',[],'YTick',[]);
%colorbar




%% AJIVE
rng(0)
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

save(['X_1_ajive.txt'],'X_1_ajive','-ascii');
save(['X_2_ajive.txt'],'X_2_ajive','-ascii');

save(['C_1_ajive.txt'],'C_1_ajive','-ascii');
save(['C_2_ajive.txt'],'C_2_ajive','-ascii');

save(['D_1_ajive.txt'],'D_1_ajive','-ascii');
save(['D_2_ajive.txt'],'D_2_ajive','-ascii');



%% JIVE: see simulation3_JIVE_and_RJIVE.R for running JIVE algorithm
X_1_jive=load('X_1_jive.txt'); 
X_2_jive=load('X_2_jive.txt');

C_1_jive=load('C_1_jive.txt');
C_2_jive=load('C_2_jive.txt');

D_1_jive=load('D_1_jive.txt');
D_2_jive=load('D_2_jive.txt');

%% R.JIVE: see simulation3_JIVE_and_RJIVE.R for running R.JIVE algorithm
X_1_rjive=load('X_1_rjive.txt'); 
X_2_rjive=load('X_2_rjive.txt');

C_1_rjive=load('C_1_rjive.txt');
C_2_rjive=load('C_2_rjive.txt');

D_1_rjive=load('D_1_rjive.txt');
D_2_rjive=load('D_2_rjive.txt');

%% DISCO-SCA
%%normalize each data matrix to unit sum of squares
norm_Y_1=norm(Y_1,'fro');
Y_1_smat=Y_1/norm_Y_1;

norm_Y_2=norm(Y_2,'fro');
Y_2_smat=Y_2/norm_Y_2;

[p_1,nn]=size(Y_1);
[p_2,~]=size(Y_2);
%By the instruction given in the article 
%"Performing DISCO-SCA to search for distinctive and common information in linked data",
%the data file concatenates all the data matrices
%in the form of (n+1)*(p_1+...+p_k) 
%with the first row gives the labels for each data matrix.
con_data=[Y_1_smat',Y_2_smat'];
row1=[ones(1,p_1),2*ones(1,p_2)];
con_data=[row1;con_data];
save('concatenated_data_c_smat.txt','con_data','-ascii');

%package obtained from http://ppw.kuleuven.be/okp/software/disco-sca/
addpath('E:\Hai Windows\work\softwares\matlab_toolbox\DISCO-SCA_Source');
%%%First to decide the number of total common and distinct components
rng(0)
DISCO_SCA
%%%We choose total number=8 by DISCOrankselection_ScreePlots.png
%%%Next, we run the DISCO-SCA
rng(0)
DISCO_SCA
%see the parameter settings in DISCOsettings.JPG

%choose 7 distinctive components according to
%DISCOresult_status_components.png
result=load('DISCOresult_7.mat')
result.DISCO_SCA.AnalysisOptions.Target
%     1     1     0     0     0     0     0     1
%     0     0     1     1     1     1     1     1
%The 8-th column is all 1, which indicates the 8-th component is common.

result.DISCO_SCA.AnalysisOutput

T=result.DISCO_SCA.AnalysisOutput.RotatedScores;%300*8
P=result.DISCO_SCA.AnalysisOutput.RotatedLoadings;%1200*8


X_1_DISCO=(T*P(1:p_1,:)')'*norm_Y_1;
X_2_DISCO=(T*P((p_1+1):(p_1+p_2),:)')'*norm_Y_2;

C_1_DISCO=(T(:,8)*P(1:p_1,8)')'*norm_Y_1;
C_2_DISCO=(T(:,8)*P((p_1+1):(p_1+p_2),8)')'*norm_Y_2;

D_1_DISCO=X_1_DISCO-C_1_DISCO;
D_2_DISCO=X_2_DISCO-C_2_DISCO;

save(['X_1_DISCO.txt'],'X_1_DISCO','-ascii');
save(['X_2_DISCO.txt'],'X_2_DISCO','-ascii');

save(['C_1_DISCO.txt'],'C_1_DISCO','-ascii');
save(['C_2_DISCO.txt'],'C_2_DISCO','-ascii');

save(['D_1_DISCO.txt'],'D_1_DISCO','-ascii');
save(['D_2_DISCO.txt'],'D_2_DISCO','-ascii');

%% COBE
rng(0)
addpath('E:\Hai Windows\work\softwares\matlab_toolbox\demo_CIFA');
%denoise 

[U1 D1 V1]=svd(Y_1_smat,'econ');
[U2 D2 V2]=svd(Y_2_smat,'econ');

Y_1_smat_denoise=U1(:,1:3)*D1(1:3,1:3)*V1(:,1:3)';
Y_2_smat_denoise=U2(:,1:5)*D2(1:5,1:5)*V2(:,1:5)';

Y{1}=Y_1_smat_denoise';
Y{2}=Y_2_smat_denoise';

Ac=cobe(Y);%No common basis found

C_1_cobe=zeros(p_1,nn);
C_2_cobe=zeros(p_2,nn);

%{
C1_cobe=X1_noisy_c_smat_denoise*Ac*Ac'*norm_X1_noisy_c;
C2_cobe=X2_noisy_c_smat_denoise*Ac*Ac'*norm_X2_noisy_c;

X1_cobe=X1_noisy_c_smat_denoise*norm_X1_noisy_c;
X2_cobe=X2_noisy_c_smat_denoise*norm_X2_noisy_c;

D1_cobe=X1_cobe-C1_cobe;
D2_cobe=X2_cobe-C2_cobe;
%}



%% OnPLS
%see simulation3_OnPLS.py for running OnPLS algorithm
X_1_onpls=load('X_1_onpls.txt'); 
X_2_onpls=load('X_2_onpls.txt');

C_1_onpls=load('C_1_onpls.txt');
C_2_onpls=load('C_2_onpls.txt');

D_1_onpls=load('D_1_onpls.txt');
D_2_onpls=load('D_2_onpls.txt');

%% plots of other methods' common matrices
[ha, pos] =tight_subplot(2,6,[.03 .03],[.05 .02],[.05 .02]) ;
%subplot(2,5,1)
 axes(ha(1)); 
      imagesc(C_1_jive)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[]);
%colorbar

%subplot(2,5,1)
 axes(ha(2)); 
      imagesc(C_1_rjive)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar
 
 
%subplot(2,5,2)
 axes(ha(3)); 
      imagesc(C_1_ajive)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,5,3)
 axes(ha(4)); 
   imagesc(C_1_onpls)
 colormap(darkb2r(-scale,scale ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,5,4)
 axes(ha(5)); 
      imagesc(C_1_DISCO)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar
    
%subplot(2,5,5)
 axes(ha(6)); 
      imagesc(C_1_cobe)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,5,6)
 axes(ha(7)); 
      imagesc(C_2_jive)
 colormap(darkb2r(  -scale,scale   ))

%colorbar
    
%subplot(2,5,6)
 axes(ha(8)); 
      imagesc(C_2_rjive)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar

%subplot(2,5,7)
 axes(ha(9)); 
      imagesc(C_2_ajive)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar

%subplot(2,5,8)
 axes(ha(10)); 
      imagesc(C_2_onpls)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar
    


%subplot(2,5,9)
 axes(ha(11)); 
   imagesc(C_2_DISCO)
 colormap(darkb2r(-scale,scale ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar


%subplot(2,5,10)
 axes(ha(12)); 
      imagesc(C_2_cobe)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
%colorbar

%Horizontal Colorbar 
 colormap(darkb2r(  -scale,scale   ))
 colorbar
 
%%Block GDFM maps
phi_y1 = load('phi_y1_gdfm.txt')';
phi_y2 = load('phi_y2_gdfm.txt')';

psi_y1 = load('psi_y1_gdfm.txt')';
psi_y2 = load('psi_y2_gdfm.txt')';

nu_y1 = load('nu_y1_gdfm.txt')';
nu_y2 = load('nu_y2_gdfm.txt')';

xi_xy1 = load('xi_xy1_gdfm.txt')';
xi_xy2 = load('xi_xy2_gdfm.txt')';

[ha, pos] =tight_subplot(2,6,[.03 .03],[.05 .02],[.05 .02]) ;
%subplot(2,5,1)
 axes(ha(1)); 
      imagesc(phi_y1)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[]);
 
  axes(ha(2)); 
      imagesc(psi_y1)
 colormap(darkb2r(  -scale,scale   ))
  set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(3)); 
      imagesc(nu_y1)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(4)); 
      imagesc(xi_xy1)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(5)); 
 %leave it blank 
 set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(6)); 
 %leave it blank 
 set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(7)); 
      imagesc(phi_y2)
 colormap(darkb2r(  -scale,scale   ))
 
 
  axes(ha(8)); 
      imagesc(psi_y2)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(9)); 
      imagesc(nu_y2)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(10)); 
      imagesc(xi_xy2)
 colormap(darkb2r(  -scale,scale   ))
 set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(11)); 
 %leave it blank 
 set(gca,'XTick',[],'YTick',[]);
 
  axes(ha(12)); 
 %leave it blank 
 set(gca,'XTick',[],'YTick',[]);
 
 
