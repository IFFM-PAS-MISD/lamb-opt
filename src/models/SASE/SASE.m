%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%          Semi Analytical Finite Element             %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is based on "A detailed study of guided wave propagation in a
% viscoelastic multilayered anisotropic plate" paper by L. Taupin et. al.
% in 9th Anglo-French Physical Acoustic Joint conference in 2011.
% Further reference: Modeling wave propagation in damped waveguides of
% arbitrary cross-section, by I. Bartoli et. al. (2006)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs: layup, material properties
% layup angles are as given in the paper: +Z = 0 deg. and then
% counter-clockwise 
% For material properties transformation, +Y = 0 deg. and then
% counter-clockwise -- for stacking direction +X.
% Hence, layup angles have to be changed. rot_angles = layup + 90 deg.

%%
% clear
% clc

function [mode_shapes, wave_number_real, wave_number_imag] = SASE(nele_layer,np,beta,freq)

freq = freq*1e3;

%% Input
SASEinput

%% Transform material properties

C = cell(nlayers,1);

for ii=1:nlayers
    theta = rot_angles(ii);
    C_tmp = transprop(C0r+1i*C0i,stack_dir,1*theta);
    C{ii} = C_tmp;
end


%% Get K and M
% Get global stiffness and mass matrices (independent of beta and freq.)
[K, M] = get_KM(nele_layer,np,h,C,rho);

%% Loop over propagating directions (beta)

nbeta = length(beta);

wave_number_real = zeros(length(M),nbeta);
wave_number_imag = zeros(length(M),nbeta);
mode_shapes = cell(nbeta,1);

tic
for ii=1:nbeta
%     disp('Angle'); disp(beta(ii));
    [mode_shapes{ii}, wave_number_real(:,ii), wave_number_imag(:,ii)] = get_k(K,M,beta(ii),freq);
end
toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




