function [wavenumber,CG,FREQ] = main_SASE2_damp(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,C11i,C12i,C13i,C22i,C23i,C33i,C44i,C55i,C66i,layup,h,fmin,fmax,number_of_frequency_points,beta,stack_dir,np,nele_layer)
% MAIN_SASE_DAMP   Dispersion curves of Lamb wave modes 
%    includes symmetric, antisymmetric and shear horizontal modes 
%    sweep over frequencies 
%    complex values of elastic constant matrix (attenuation considered) 
% 
% Syntax: [wavenumber,CG,FREQ] = main_SASE2_damp(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,C11i,C12i,C13i,C22i,C23i,C33i,C44i,C55i,C66i,layup,h,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer) 
% 
% Inputs:
%    rho - mass density of composite layers, double, Units: kg/m^3
%    C11,C12,C13,C22,C23,C33,C44,C55,C66 - elastic constants, double, Units: Pa
%                                          dimensions [1,1] or [m,1] where
%                                          m is the number of composite layers
%    C11i,C12i,C13i,C22i,C23i,C33i,C44i,C55i,C66i - imaginary part of elastic constants, double, Units: Pa
%                                          dimensions [1,1] or [m,1] where
%                                          m is the number of composite layers
%    layup - composite layup (angles of reinforcing fibres), dimensions [1, m], Units: deg 
%    h - thickness of each layer, double, dimensions [m, 1], Units: m
%    wavenumber_min - minimum value of wavenumber for dispersion curves,double, dimensions [number_of_angles,1], Units: rad/m
%    wavenumber_max - maximum value of wavenumber for dispersion curves,double, dimensions [number_of_angles,1], Units: rad/m
%    number_of_wavenumber_points - number of wavenumber points, integer
%    beta - vector of wave propagation angles under analysis, double,dimensions [1,number_of_angles]
%    stack_dir - stacking direction (1,2 or 3)
%    np - order of elements (3<=np<=5)
%    nele_layer - number of spectral elements per ply, integer, default=1
% 
% Outputs: 
%    wavenumber - vector of wavenumbers, dimensions [1, number_of_wavenumber_points], Units: rad/m 
%    CG - matrix of group velocities, double, dimensions [number_of_modes,number_of_wavenumber_points,number_of_angles], Units: m/s 
%    FREQ - matrix of frequencies, double, dimensions [number_of_modes,number_of_wavenumber_points,number_of_angles], Units: Hz 
% 
% Example: 
%    [wavenumber,CG,FREQ] = main_SASE2_damp(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,C11i,C12i,C13i,C22i,C23i,C33i,C44i,C55i,C66i,layup,h,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer) 
%    [wavenumber,CG,FREQ] = main_SASE2_damp(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,C11i,C12i,C13i,C22i,C23i,C33i,C44i,C55i,C66i,[0,90,0,90],h,wavenumber_max,512,[0,30,60,90],1,3,1) 
% 
% Other m-files required: get_om 
% Subfunctions: none 
% MAT-files required: none 
% See also: 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is based on "A detailed study of guided wave propagation in a
% viscoelastic multilayered anisotropic plate" paper by L. Taupin et. al.
% in 9th Anglo-French Physical Acoustic Joint conference in 2011.
% Further reference: Modeling wave propagation in damped waveguides of
% arbitrary cross-section, by I. Bartoli et. al. (2006)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% layup angles are as given in the paper: +Z = 0 deg. and then
% counter-clockwise 
% For material properties transformation, +Y = 0 deg. and then
% counter-clockwise -- for stacking direction +X.
% Hence, layup angles have to be changed. rot_angles = layup - 90 deg.

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

nlayers = length(layup);
rot_angles = layup - 90; % set x axis horizontal 
number_of_angles = length(beta);
frequency_step=(fmax-fmin)/(number_of_frequency_points-1); % frequency step [Hz]

%% Transform material properties
C = cell(nlayers,1);
C0r = zeros(6,6);
C0i = zeros(6,6);
for ii=1:nlayers
    theta = rot_angles(ii);
    if(length(C11) == 1)
        C0r(1,1) = C11; C0r(1,2)=C12;C0r(2,1)=C12; C0r(1,3)=C13;C0r(3,1)=C13;C0r(2,2)=C22;C0r(2,3)=C23;C0r(3,2)=C23;C0r(3,3)=C33;C0r(4,4)=C44;C0r(5,5)=C55;C0r(6,6)=C66;
        C0i(1,1) = C11i; C0i(1,2)=C12i;C0i(2,1)=C12i; C0i(1,3)=C13i;C0i(3,1)=C13i;C0i(2,2)=C22i;C0i(2,3)=C23i;C0i(3,2)=C23i;C0i(3,3)=C33i;C0i(4,4)=C44i;C0i(5,5)=C55i;C0i(6,6)=C66i;
    else
        C0r(1,1) = C11(ii); C0r(1,2)=C12(ii); C0r(2,1)=C12(ii);C0r(1,3)=C13(ii);C0r(3,1)=C13(ii);C0r(2,2)=C22(ii);C0r(2,3)=C23(ii);C0r(3,2)=C23(ii);C0r(3,3)=C33(ii);C0r(4,4)=C44(ii);C0r(5,5)=C55(ii);C0r(6,6)=C66(ii);
        C0i(1,1) = C11i(ii); C0i(1,2)=C12i(ii); C0i(2,1)=C12i(ii);C0i(1,3)=C13i(ii);C0i(3,1)=C13i(ii);C0i(2,2)=C22i(ii);C0i(2,3)=C23i(ii);C0i(3,2)=C23i(ii);C0i(3,3)=C33i(ii);C0ri(4,4)=C44i(ii);C0i(5,5)=C55i(ii);C0r(6,6)=C66(ii);
    end
    %C_tmp = transform_prop_3D(C0r,stack_dir,theta); % no damping, real elastic constants
    %C_tmp = transform_prop_3D(C0r+1i*C0i,stack_dir,theta); % damping, complex elastic constants
    C_tmp = transform_prop_3D(C0r-1i*C0i,stack_dir,theta); % damping, complex elastic constants
    C{ii} = C_tmp;
end

%% Get K and M
% Get global stiffness and mass matrices (independent of beta and wavenumber)
[K, M] = get_KM(nele_layer,np,h,C,rho);
%% dispersion curves
number_of_modes=length(M);
FREQ = zeros(number_of_frequency_points,number_of_angles);
CG = zeros(number_of_modes,number_of_frequency_points,number_of_angles);
wavenumber = zeros(number_of_modes,number_of_frequency_points,number_of_angles);
for k=1:number_of_frequency_points
    FREQ(k,:) = fmin+(k-1)*frequency_step;
end
%% loop over angles
for j=1:length(beta)
    %fprintf('SASE dispersion curves at angle: %2.1f\n', beta(j));
    wave_number_real = zeros(number_of_modes,number_of_frequency_points);
    %om_imag = zeros(number_of_modes,number_of_frequency_points);
    %cg = zeros(number_of_modes,number_of_frequency_points);
    for k=1:number_of_frequency_points
        %[cg(:,k),~, om_real(:,k), ~] = get_om(K,M,beta(j),wavenumber(k,j));    
        [mode_shapes, wave_number_real(:,k), wave_number_imag] = get_k(K,M,beta(j),FREQ(k,j));
    end
    %% mode-tracing
    % Taylor approximation method
    %[cg_new,om_new] = mode_tracing(cg,om_real,wavenumber_step(j));
    %[cg_new,om_new] = mode_tracing_pade(cg,om_real,wavenumber_step(j));
      %cg_new=cg;
  
    wavenumber(:,:,j) = wave_number_real;
    %CG(:,:,j) = cg_new;
end
%---------------------- END OF CODE---------------------- 

% ================ [main_SASE2_damp.m] ================  
