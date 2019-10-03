%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                       get_om.m                       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves generalized eigenvalue problem 
% [A-om^2*M] *eig_vect = 0
% where in this case eigen values are om^2
% om are angular frequency components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cg,QR,om_real,om_imag] = get_om(K,M,beta,wavenumber)

b = beta*pi/180;

K11 = K{1,1};   K12 = K{1,2};   K13 = K{1,3};
K21 = K{2,1};   K22 = K{2,2};   K23 = K{2,3};
K31 = K{3,1};   K32 = K{3,2};   K33 = K{3,3};

cb = cos(b);
sb = sin(b);
% cb=1;
% sb=0;
m = length(M);

T = sparse(m,1);
T(1:3:end) = 1i;
T(2:3:end) = 1;
T(3:3:end) = 1;
% T(1:3:end) = 1;
% T(2:3:end) = 1i;
% T(3:3:end) = 1i;

T=diag(T);

A = wavenumber^2*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32)+1i*wavenumber*T'*(-cb*K13 - sb*K21 + sb*K12 + cb*K31)*T + K11;

%% Standard Eigenvalue problem

[eig_vect, eig_val,eig_vect_r] = eig(A,M);
%[eig_vect, eig_val,eig_vect_r] = eig(K11,M); % cut-off frequencies
eig_val = (diag(eig_val));


[om_real, IX_real] = sort(abs(real(eig_val(1:1:end))),'ascend');
om_real = sqrt(om_real);

%[om_imag, IX_imag] = sort(abs(imag(eig_val(1:1:end))),'ascend');
[om_imag] = sqrt(abs(imag(eig_val(1:1:end))));
% om_real
% om_imag

QR = (eig_vect(:,IX_real));
QL = (eig_vect_r(:,IX_real))';
%QL = (eig_vect_r(IX_real,:))';

 for k=1:length(M)
     if(~isreal(eig_val(IX_real(k))))
         QR(:,k) = abs(QR(:,k));
         QL(k,:) = abs(QL(k,:));
%          QR(:,k) = real(QR(:,k));
%          QL(k,:) = real(QL(k,:));
%          QR(:,k) = real(conj(QR(:,k)));
%          QL(k,:) = real(conj(QL(k,:)));
     end
 end
%     else
%         QR(:,k) = (QR(:,k));
%     end
% end
%QR = (eig_vect(:,IX_real));
%QR = real(eig_vect(:,IX_real));

%%
% group velocity computation
Aprim = 2*wavenumber*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32) +1i*T'*(-cb*K13 - sb*K21 + sb*K12 + cb*K31)*T;

% if (~isreal(QR)) 
%     disp('complex'); 
%     return 
% end
cg = zeros(length(M),1);
for k=1:length(M)
    QL=QR(:,k)';
    cg(k,1) =  (QL*Aprim*QR(:,k))./(2*om_real(k)*QL*M*QR(:,k));
    %cg(k,1) =  (QL(k,:)*Aprim*QR(:,k))./(2*om_real(k)*QL(k,:)*M*QR(:,k));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


