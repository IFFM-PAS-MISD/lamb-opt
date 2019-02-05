function [C11,C12,C13,C22,C23,C33,C44,C55,C66] = lamina_3D(e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23)
% LAMINA_3D   Elastic constants of composite lamina in terms of principal material directions
%    Elastic constants of a fibre reinforced composite material  
%    3D elasticity 
% 
% Syntax: [C11,C12,C13,C22,C23,C33,C44,C55,C66] = lamina_3D(e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23)
% 
% Inputs: 
%    rho - density of composite lamina, double, Units: kg/m^3 
%    e11 - Young modulus of composite lamina parallel to fibres, double, Units: Pa 
%    e22 - Young modulus of composite lamina perpendicular to fibres (in-plane), double, Units: Pa 
%    e33 - Young modulus of composite lamina perpendicular to fibres (out-of-plane), double, Units: Pa 
%    ni12, ni13, ni21, ni23, ni31, ni32 - Poisson ratios, double, Units: -
%    g12, g13, g23 - shear modulus, double, Units: Pa
%
%    inputs takes also vectors of dimensions [m,1]
%    where m is the number of layers in composite laminate
% 
% Outputs: 
%    C11,C12,C13,C22,C23,C33,C44,C55,C66 - elastic constants, double, Units: Pa 
%
%    {sigma_1} = |C11 C12 C13  0   0   0  | {epsilon_1}
%    {sigma_2}   |C12 C22 C23  0   0   0  | {epsilon_2}
%    {sigma_3}   |C13 C23 C33  0   0   0  | {epsilon_3}
%    {sigma_4}   | 0   0   0  C44  0   0  | {gamma_23}
%    {sigma_5}   | 0   0   0   0  C55  0  | {gamma_31}
%    {sigma_6}   | 0   0   0   0   0  C66 | {gamma_12}
% 
% Example: 
%    [C11,C12,C13,C22,C23,C33,C44,C55,C66] = lamina_3D(e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23)
% 
% Other m-files reCuired: none 
% Subfunctions: none 
% MAT-files reCuired: none 
% See also: HOMOGENIZATION 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sierakowski page 46 eC. 2.33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delta=(1-ni12.*ni21-ni13.*ni31-ni12.*ni23.*ni31-ni13.*ni21.*ni32-ni23.*ni32);
C11=e11.*(1-ni23.*ni32)./Delta;
C22=e22.*(1-ni31.*ni13)./Delta;
C33=e33.*(1-ni12.*ni21)./Delta;
C44=g23;
C55=g13;
C66=g12;
C12=(ni21+ ni31.*ni23).*e11./Delta;
C13=(ni31+ ni21.*ni32).*e11./Delta;
C23=(ni32+ ni12.*ni31).*e22./Delta;

%---------------------- END OF CODE---------------------- 

% ================ [lamina_3D.m] ================  
