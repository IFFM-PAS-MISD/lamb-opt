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

function [cg,mode_shapes, om_real, om_imag] = SASE_om(K,M,beta,wavenumber)



%% Loop over propagating directions (beta)

number_of_angles = length(beta);

om_real = zeros(length(M),number_of_angles);
om_imag = zeros(length(M),number_of_angles);
cg = zeros(length(M),number_of_angles);
mode_shapes = cell(number_of_angles,1);


for ii=1:number_of_angles
    [cg(:,ii), mode_shapes{ii}, om_real(:,ii), om_imag(:,ii)] = get_om(K,M,beta(ii),wavenumber);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




