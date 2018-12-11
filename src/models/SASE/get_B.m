%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                       get_B.m                       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gets B1, B2 and B3 matrices for a given element number
% of order (np) and at given integration point (ii_th integration point)
% Given: Jacobian at ii_th integration point (J) and matrices of shape
% functions (SF) and their derivatives (SFD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B, N] = get_B(ii,np,J,SF,SFD)

intp = ii;
nnodes = np+1;

B = cell(3,1);

% Declare and populate Lx, Ly and Lz -- ndof = 3 (do NOT change)
Lx = zeros(6,3);
Ly = zeros(6,3);
Lz = zeros(6,3);

Lx(1,1) = 1;    Lx(5,3) = 1;    Lx(6,2) = 1;
Ly(2,2) = 1;    Ly(4,3) = 1;    Ly(6,1) = 1;
Lz(3,3) = 1;    Lz(4,2) = 1;    Lz(5,1) = 1;

% Values and derivatives of shape functions in parent domain

N = zeros(3,nnodes*3);
N_x = zeros(3,nnodes*3);

for inode=1:nnodes
    t = (inode-1)*3;
    
    N(1,t+1) = SF(inode,intp);
    N(2,t+2) = SF(inode,intp);
    N(3,t+3) = SF(inode,intp);
    
    N_x(1,t+1) = SFD(inode,intp);
    N_x(2,t+2) = SFD(inode,intp);
    N_x(3,t+3) = SFD(inode,intp);
    
end

if(abs(J)<=1e-10)
    disp('J ~ 0')
end

N_x = N_x / J;
            
B1 = Lx * N_x;
B2 = Ly * N;
B3 = Lz * N;

B{1} = B1;
B{2} = B2;
B{3} = B3;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N = [psi1      0       0       psi2        0       0       psi3        0       0;
%        0       psi1    0          0        psi2    0          0        psi3    0;
%        0       0       psi1       0        0       psi2       0        0       psi3];

% N_x = 1/J*[psi1_zi      0           0           psi2_zi         0           0           psi3_zi     0           0;
%                 0       psi1_zi     0           0               psi2_zi     0           0           psi3_zi     0;
%                 0       0           psi1_zi     0               0           psi2_zi   	0           0           psi3_zi];