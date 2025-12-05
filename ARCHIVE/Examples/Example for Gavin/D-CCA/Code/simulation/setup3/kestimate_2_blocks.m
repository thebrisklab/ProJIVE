%adapted from function kestimate_3_blocks_SS in http://www.barigozzi.eu/BHS_final_codes.zip
%
function [chi_xy1, chi_y1, phi_y1, psi_y1, nu_y1, xi_xy1,   chi_xy2, chi_y2, phi_y2, psi_y2, nu_y2, xi_xy2] = kestimate_2_blocks(y1,y2,q12,q1,q2,M)  


x = [y1,y2];
[T, Nx] = size(x);


[T, Ny1] = size(y1);
[T, Ny2] = size(y2);


w = M;
W = 2*w+1;

% compute covariances x
Sx = zeros(W,Nx*Nx);
for k = 1:w+1,
     Sk = (x(k:T,:))'*(x(1:T+1-k,:))/(T-k);
     S_k = Sk';
     Sx(w-k+2,:) = S_k(:)'; 
     Sx(w+k,:) = Sk(:)';     
end
% compute the spectral matrix in W points (S)
Factor = exp(-sqrt(-1)*(-w:w)'*(0:2*pi/W:4*pi*w/W)); %  (2w+1) x  (2w+1) 
Sx = diag(triang(W))*Sx(1:W,:);
Sxx = reshape(Sx'*Factor,Nx,Nx*W);

% compute the first nfactors eigenvalues and 
% eigenvectors  for all points (D, A) and the K function
Pa1 = zeros(Ny1,Nx,W);
Pa2 = zeros(Ny2,Nx,W);


Sx = Sxx(:,1:Nx);

% projection furter 1st block  on q principal components of 3 blocks

i1 = [1:Ny1];
i2 = [Ny1+1:Ny1+Ny2];


for j = 1:w+1,

    Sx = Sxx(:,(j-1)*Nx+1:j*Nx);
    [Ax,Dx] = eigs(Sx, q12);

    Dx = diag(diag(real(Dx)));
    Dxi = diag(diag(real(Dx)).^(-1));
    
    % 1st projection
    Pa1(:,:,j) = Ax(i1,:)*Ax';
    Pa2(:,:,j) = Ax(i2,:)*Ax';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
    [A121,D121] = eigs(Ax(i1,:)*Dx*Ax(i1,:)', q1);
    [A122,D122] = eigs(Ax(i2,:)*Dx*Ax(i2,:)', q2);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pc121(:,:,j) =  A121*A121'*Ax(i1,:)*Ax';
    Pc122(:,:,j) =  A122*A122'*Ax(i2,:)*Ax';


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  12
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    D121i2 = diag(diag(real(D121)).^(-0.5));
    D122i = diag(diag(real(D122)).^(-1));
    
    SV12 = D121i2*(A121'*Ax(i1,:)*Dx*(  A122'*Ax(i2,:)  )')*...
         D122i*...
         (D121i2*(A121'*Ax(i1,:)*Dx*(  A122'*Ax(i2,:)  )'))';
     
     if q1+q2-q12
         
         [A5,D5] = eigs(SV12, q1+q2-q12);
         
         %if abs(sum(sum(A5'*A5-eye(q1+q2-q12)))) > 1e-05             
         %    A5 = eye(q1+q2-q12);
         %end
         
         Base = A5'*D121i2*(A121'*Ax(i1,:));
         Base12 = Base*Ax';
         
         
         Pd1w12(:,:,j) = A121*A121'*Ax(i1,:)*Dx*...
             (A121'*Ax(i1,:))'*D121i2'*A5*...
             Base12;
         
         
         Pd2w12(:,:,j) = A122*A122'*Ax(i2,:)*Dx*...
             (A121'*Ax(i1,:))'*D121i2'*A5*...
             Base12;
         
     end
if j>1

    Pa1(:,:,W+2-j) = conj(Pa1(:,:,j));
    Pa2(:,:,W+2-j) = conj(Pa2(:,:,j));

    Pc121(:,:,W+2-j) = conj(Pc121(:,:,j));
    Pc122(:,:,W+2-j) = conj(Pc122(:,:,j));
 
    
    if q1+q2-q12      
        Pd1w12(:,:,W+2-j) = conj(Pd1w12(:,:,j));
        Pd2w12(:,:,W+2-j) = conj(Pd2w12(:,:,j));    
    end
    
      
    
end 

end



chi_xy1 = zeros(T+2*w,Ny1);
chi_xy2 = zeros(T+2*w,Ny2);

chi_y1 = zeros(T+2*w,Ny1);
chi_y2 = zeros(T+2*w,Ny2);

phi_y1 = zeros(T+2*w,Ny1);
phi_y2 = zeros(T+2*w,Ny2);

for j = 1:Nx
  U_chi_xy1 = real(conj(Factor)*squeeze(Pa1(:,j,:)).'/W);
  chi_xy1  = chi_xy1 + conv2(U_chi_xy1,x(:,j));
    
  U_chi_xy2 = real(conj(Factor)*squeeze(Pa2(:,j,:)).'/W);
  chi_xy2  = chi_xy2 + conv2(U_chi_xy2,x(:,j));
  
  
  U_chi_y1 = real(conj(Factor)*squeeze(Pc121(:,j,:)).'/W);
  chi_y1  = chi_y1 + conv2(U_chi_y1,x(:,j));
                   
  U_chi_y2 = real(conj(Factor)*squeeze(Pc122(:,j,:)).'/W);
  chi_y2  = chi_y2 + conv2(U_chi_y2,x(:,j));
  
  if q1+q2-q12
    U_phi_y1 = real(conj(Factor)*squeeze(Pd1w12(:,j,:)).'/W);
    phi_y1  = phi_y1 + conv2(U_phi_y1,x(:,j));
    
    U_phi_y2 = real(conj(Factor)*squeeze(Pd2w12(:,j,:)).'/W);
    phi_y2  = phi_y2 + conv2(U_phi_y2,x(:,j));
  end
end

chi_xy1 = chi_xy1(w + 1:T + w ,:);
chi_xy2 = chi_xy2(w + 1:T + w ,:);

chi_y1 = chi_y1(w + 1:T + w ,:);
chi_y2 = chi_y2(w + 1:T + w ,:);

phi_y1 = phi_y1(w + 1:T + w ,:);
phi_y2 = phi_y2(w + 1:T + w ,:);

psi_y1 = chi_y1 - phi_y1;
psi_y2 = chi_y2 - phi_y2;


nu_y1 = chi_xy1 - chi_y1;
nu_y2 = chi_xy2 - chi_y2;

xi_xy1 = y1 - chi_xy1;
xi_xy2 = y2 - chi_xy2;





