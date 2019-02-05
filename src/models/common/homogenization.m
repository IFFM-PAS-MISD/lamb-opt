function [rho,e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23] = homogenization(rhom,rhof,em,ef,nim,nif,vol)
% HOMOGENIZATION   Homogenization of fibre reinforced composite properties 
%    homogenize material properties of composite with unidirectional fibres 
%    based on its constituents: matrix (m) and fibres (f) 
% 
% Syntax: [rho,e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23] = homogenization(rhom,rhof,em,ef,nim,nif,vol) 
% 
% Inputs: 
%    rhom - density of matrix, double, Units: kg/m^3 
%    rhof - density of fibres, double, Units: kg/m^3  
%    em - Young modulus of matrix, double, Units: Pa 
%    ef - Young modulus of fibres, double, Units: Pa 
%    vol - volume fraction of reinforcing fibres in range [0:1], double, Units: - 
%
%    inputs takes also vectors of dimensions [m,1]
%    where m is the number of layers in composite laminate
% 
% Outputs: 
%    rho - density of composite lamina, double, Units: kg/m^3 
%    e11 - Young modulus of composite lamina parallel to fibres, double, Units: Pa 
%    e22 - Young modulus of composite lamina perpendicular to fibres (in-plane), double, Units: Pa 
%    e33 - Young modulus of composite lamina perpendicular to fibres (out-of-plane), double, Units: Pa 
%    ni12, ni13, ni21, ni23, ni31, ni32 - Poisson ratios, double, Units: -
%    g12, g13, g23 - shear modulus, double, Units: Pa
% 
% Example: 
%    [rho,e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23] = homogenization(rhom,rhof,em,ef,nim,nif,vol)
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: LAMINA_3D 
% 
% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

rho=rhof.*vol+rhom.*(1-vol);
gm=em./(1+nim)/2; gf=ef./(1+nif)/2;
e11=ef.*vol+em.*(1-vol);
e22=((ef+em)+(ef-em).*vol)./((ef+em)-(ef-em).*vol).*em;
e33=e22;
ni12=nif.*vol+nim.*(1-vol);
ni13=ni12;
ni21=e22./e11.*ni12;
ni31=e33./e11.*ni13;
ni23=nif.*vol+nim.*(1-vol) .*(1.+nim-ni12.*em./e11) ./(1.-nim.^2 +nim.*ni12.*em./e11);
ni32=e33./e22.*ni23;
g12=((gf+gm)+(gf-gm).*vol)./((gf+gm)-(gf-gm).*vol).*gm;
g23=e22./2./(1+ni23);
g13=g12;

%---------------------- END OF CODE---------------------- 

% ================ [homogenization.m] ================  
