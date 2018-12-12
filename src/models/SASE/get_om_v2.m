%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                       get_om.m                       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves generalized eigenvalue problem 
% [A-eig_val^2] *eig_vect = 0
% where in this case eig_val are frequency components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cg,QR,om_real,om_imag] = get_om_v2(K,M,wavenumber)


K1 = K{1,1};   K2 = K{1,2} -  K{2,1};   K3 = K{2,2};  


m = length(M);

T = sparse(m,1);
T(1:3:end) = 1i;
T(2:3:end) = 1;
T(3:3:end) = 1;

A = K1 +1i *wavenumber* (T'*K2*T) + wavenumber^2*K3;

%% Generalized Eigenvalue problem

[eig_vect, eig_val] = eig(A,M);
eig_val = sqrt(diag(eig_val));

[om_real, IX_real] = sort((real(eig_val(1:1:end))),'ascend');

[om_imag, IX_imag] = sort((imag(eig_val(1:1:end))),'ascend');


QR = (eig_vect(:,IX_real));

%%
% group velocity computation
Aprim = 2*wavenumber*K3 + 1i * (T'*K2*T);

cg = zeros(length(M),1);
for k=1:length(M)
    QL=QR(:,k)';
    cg(k,1) =  (QL*Aprim*QR(:,k))./(2*om_real(k)*QL*M*QR(:,k));   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


