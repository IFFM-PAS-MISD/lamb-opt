%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                       get_k.m                       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves generalized eigenvalue problem 
% A*eig_vect = B*eig_vect*eig_val
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mode_shapes, wave_number_real, wave_number_imag] = get_k(K,M,beta,freq)

b = beta*pi/180;
omega = 2*pi*freq;

K11 = K{1,1};   K12 = K{1,2};   K13 = K{1,3};
K21 = K{2,1};   K22 = K{2,2};   K23 = K{2,3};
K31 = K{3,1};   K32 = K{3,2};   K33 = K{3,3};

cb = cos(b);
sb = sin(b);


A = [zeros(size(M))     K11-omega^2*M;
    K11-omega^2*M       -1i*(cb*K13 + sb*K21 - sb*K12 - cb*K31)];

D = [K11-omega^2*M      zeros(size(M));
    zeros(size(M))      -(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32)];


%% Generalized Eigenvalue problem

[eig_vect, eig_val] = eig(A,D);
eig_val = sort(diag(eig_val));


%[wave_number_real, IX_real] = sort(abs(real(eig_val(1:2:end))));

[wave_number_imag, IX_imag] = sort(abs(imag(eig_val(1:2:end))));
%[wave_number_imag] = abs(imag(eig_val(1:2:end)));

[wave_number_real] = abs(real(eig_val(1:2:end)));
wave_number_real=wave_number_real(IX_imag);
mode_shapes = eig_vect(:,IX_imag);
%mode_shapes = eig_vect;

% eig_val = sort(diag(eig_val));
% 
% wave_number_real = sort(abs(real(eig_val(1:2:end))));
% wave_number_imag = sort(abs(imag(eig_val(1:2:end))));
%%
% group velocity computation
%Aprim1 = 2*wave_number_real(1)*(sb^2*K22 + cb^2*K33 - cb*sb*K23 - cb*sb*K32) -1i *(cb*K13 + sb*K21 - sb*K12 - cb*K31);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


