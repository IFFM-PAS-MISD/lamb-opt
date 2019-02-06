function [wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer)
% MAIN_SASE   Dispersion curves of Lamb wave modes 
%    includes symmetric, antisymmetric and shear horizontal modes 
%    sweep over wavenumbers 
%    only real values of elastic constant matrix (no attenuation) 
% 
% Syntax: [wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer) 
% 
% Inputs:
%    rho - mass density of composite layers, double, Units: kg/m^3
%    C11,C12,C13,C22,C23,C33,C44,C55,C66 - elastic constants, double, Units: Pa
%                                          dimensions [1,1] or [m,1] where
%                                          m is the number of composite layers
%    layup - composite layup (angles of reinforcing fibres), dimensions [1, m], Units: deg 
%    h - thickness of each layer, double, dimensions [m, 1], Units: m
%    wavenumber_min - minimum value of wavenumber for dispersion curves,double, dimensions [number_of_angles,1], Units: 1/m
%    wavenumber_max - maximum value of wavenumber for dispersion curves,double, dimensions [number_of_angles,1], Units: 1/m
%    number_of_wavenumber_points - number of wavenumber points, integer
%    beta - vector of wave propagation angles under analysis, double,dimensions [1,number_of_angles]
%    stack_dir - stacking direction (1,2 or 3)
%    np - order of elements (3<=np<=5)
%    nele_layer - number of spectral elements per ply, integer, default=1
% 
% Outputs: 
%    wavenumber - vector of wavenumbers, dimensions [1, number_of_wavenumber_points], Units: 1/m 
%    CG - matrix of group velocities, double, dimensions [num_of_modes,number_of_wavenumber_points,number_of_angles], Units: m/s 
%    FREQ - matrix of frequencies, double, dimensions [num_of_modes,number_of_wavenumber_points,number_of_angles], Units: Hz 
% 
% Example: 
%    [wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer) 
%    [wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,[0,90,0,90],h,wavenumber_max,512,[0,30,60,90],1,3,1) 
% 
% Other m-files required: SASE_om 
% Subfunctions: none 
% MAT-files required: none 
% See also: 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

nlayers = length(layup);
rot_angles = layup - 90; % set x axis horizontal 
number_of_angles = length(beta);
wavenumber_step=zeros(number_of_angles,1);
for j=1:number_of_angles
    wavenumber_step(j)=(wavenumber_max(j)-wavenumber_min(j))/(number_of_wavenumber_points-1); % wavenumber step [1/m]
end
%% Transform material properties
C = cell(nlayers,1);
C0r = zeros(6,6);
for ii=1:nlayers
    theta = rot_angles(ii);
    if(length(C11) == 1)
        C0r(1,1) = C11; C0r(1,2)=C12; C0r(1,3)=C13;C0r(2,2)=C22;C0r(2,3)=C23;C0r(3,3)=C33;C0r(4,4)=C44;C0r(5,5)=C55;C0r(6,6)=C66;
    else
        C0r(1,1) = C11(ii); C0r(1,2)=C12(ii); C0r(1,3)=C13(ii);C0r(2,2)=C22(ii);C0r(2,3)=C23(ii);C0r(3,3)=C33(ii);C0r(4,4)=C44(ii);C0r(5,5)=C55(ii);C0r(6,6)=C66(ii);
    end
    C_tmp = transprop(C0r,stack_dir,theta); % no damping, real elastic constants
    C{ii} = C_tmp;
end


%% Get K and M
% Get global stiffness and mass matrices (independent of beta and wavenumber)
[K, M] = get_KM(nele_layer,np,h,C,rho);
%% dispersion curves
num_of_modes=length(M);
FREQ = zeros(num_of_modes,number_of_wavenumber_points,number_of_angles);
CG = zeros(num_of_modes,number_of_wavenumber_points,number_of_angles);
wavenumber = zeros(number_of_wavenumber_points,number_of_angles);
for k=1:number_of_wavenumber_points
    wavenumber(k,:) = wavenumber_min+(k-1)*wavenumber_step;
end
%% loop over angles
for j=1:length(beta)
    fprintf('SASE dispersion curves at angle: %2.1f\n', beta(j));
    om_real = zeros(num_of_modes,number_of_wavenumber_points);
    om_imag = zeros(num_of_modes,number_of_wavenumber_points);
    cg = zeros(num_of_modes,number_of_wavenumber_points);
    for k=1:number_of_wavenumber_points
        [cg(:,k),mode_shapes, om_real(:,k), om_imag(:,k)] = SASE_om(K,M,beta(j),wavenumber(k,j));   
    end
    %% mode-tracing
    % Taylor approximation method
    [cg_new,om_new] = mode_tracing(cg,om_real,wavenumber_step(j));
    FREQ(:,:,j) = om_new/2/pi;
    CG(:,:,j) = cg_new;
end
%---------------------- END OF CODE---------------------- 

% ================ [main_SASE.m] ================  
