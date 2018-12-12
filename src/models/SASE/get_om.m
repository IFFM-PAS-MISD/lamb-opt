%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                       get_om.m                       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves generalized eigenvalue problem 
% [A-eig_val^2] *eig_vect = 0
% where in this case eig_val are frequency components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cg,QR,om_real,om_imag] = get_om(K,M,beta,wavenumber)

b = beta*pi/180;


K11 = K{1,1};   K12 = K{1,2};   K13 = K{1,3};
K21 = K{2,1};   K22 = K{2,2};   K23 = K{2,3};
K31 = K{3,1};   K32 = K{3,2};   K33 = K{3,3};

cb = cos(b);
sb = sin(b);
m = length(M);

T = sparse(m,1);
T(1:3:end) = 1i;
T(2:3:end) = 1;
T(3:3:end) = 1;

T=diag(T);

%A = wavenumber^2*(sb^2*T'*K22*T + cb^2*T'*K33*T - cb*sb*T'*K23*T - cb*sb*T'*K32*T) - 1i*wavenumber*(cb*T'*K13*T + sb*T'*K21*T - sb*T'*K12*T - cb*T'*K31*T) + T'*K11*T;
%A = wavenumber^2*(sb^2*T'*K22*T + cb^2*T'*K33*T - cb*sb*T'*K23*T - cb*sb*T'*K32*T) + wavenumber*(cb*T'*K13*T + sb*T'*K21*T - sb*T'*K12*T - cb*T'*K31*T) + T'*K11*T;
%
%A = wavenumber^2*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32)-1i*wavenumber*(cb*T'*K13*T + sb*T'*K21*T - sb*T'*K12*T - cb*T'*K31*T) + K11;
%A = wavenumber^2*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32)+wavenumber*(cb*T'*K13*T + sb*T'*K21*T - sb*T'*K12*T - cb*T'*K31*T) + K11;

%A = wavenumber^2*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32)-1i*wavenumber*(cb*K13 + sb*K21 - sb*K12 - cb*K31) + K11;
%A = wavenumber^2*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32)-wavenumber*(cb*K13 + sb*K21 - sb*K12 - cb*K31) + K11;
A = wavenumber^2*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32)+1i*wavenumber*T'*(-cb*K13 - sb*K21 + sb*K12 + cb*K31)*T + K11;

%% Generalized Eigenvalue problem

[eig_vect, eig_val] = eig(A,M);

eig_val = sqrt(diag(eig_val));


[om_real, IX_real] = sort((real(eig_val(1:1:end))),'ascend');

[om_imag, IX_imag] = sort((imag(eig_val(1:1:end))),'ascend');


QR = (eig_vect(:,IX_real));


%%
% group velocity computation
Aprim = 2*wavenumber*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32) +1i*T'*(-cb*K13 - sb*K21 + sb*K12 + cb*K31)*T;

cg = zeros(length(M),1);
for k=1:length(M)
    QL=QR(:,k)';
    cg(k,1) =  (QL*Aprim*QR(:,k))./(2*om_real(k)*QL*M*QR(:,k));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


