function [q11, q12, q13, q22, q23, q33, q44, q55, q66] = ...
    functionalHahn(e11_m, e11_f, e22_f, ni12_m, ni12_f, ni23_f, vol)
% FUNCTIONALHAHN   homogenization of composite elastic properties 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = functionalHahn(input1,input2,input3) 
% 
% Inputs: 
%    e11_m - matrix elastic modulus E11 , dimensions [1, 1], Units: GPa
%    ni12_m - matrix Poisson's ratio , dimensions [1, 1], Units: -
%    e11_f - fiber elastic modulus E11 , dimensions [1, 1], Units: GPa
%    e22_f - fiber elastic modulus E22 , dimensions [1, 1], Units: GPa
%    ni12_f - fiber Poisson's ratio ni12, dimensions [1, 1], Units: -
%    ni23_f - fiber Poisson's ratio ni23, dimensions [1, 1], Units: -
%    vol - fiber volume fraction , dimensions [1, n], Units: -
% Outputs: 
%    q11,...q66 - homogenized mechanical properties, integer, dimensions [1, n], Units: GPa 
%    
% 
% Example: 
%    [output1,output2] = functionalHahn(input1,input2,input3) 
%    [output1,output2] = functionalHahn(input1,input2) 
%    [output1] = functionalHahn(input1,input2,input3) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Piotr Fiborek, D.Sc., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pfiborek@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 
g12_m = e11_m./(2*(1+ni12_m)); 
g12_f = e11_f./(2*(1+ni12_f)); g23_f = e22_f/(2*(1+ni23_f));
ni23_m=ni12_m;

P_f = e11_f; P_m = e11_m; eta = 1;
P = (P_f.*vol+eta.*P_m.*(1-vol))./(vol+eta.*(1-vol));
e11 = P;
   
P_f = ni12_f; P_m = ni12_m; eta = 1;
P = (P_f.*vol+eta.*P_m.*(1-vol))./(vol+eta.*(1-vol));
ni12 = P; ni13 = ni12;

P_f = 1./g12_f; P_m = 1./g12_m; eta = (1+g12_m./g12_f)/2;
P = (P_f.*vol+eta.*P_m.*(1-vol))./(vol+eta.*(1-vol));
g12 = 1./P;  g13 = g12;
   
eta = (3-4.*ni12_m+g12_m./g12_f)./(4*(1-ni12_m));
P_f = 1/g23_f;
P = (P_f.*vol+eta.*P_m.*(1-vol))./(vol+eta.*(1-vol));
g23 = 1./P;
    
K_f = e11_f./(2*(1-ni12_f)); K_m = e11_m./(2*(1-ni12_m));
P_f = 1./K_f; P_m = 1./K_m; eta = (1+g12_m./K_f)./(2*(1-ni12_m));
P = (P_f.*vol+eta.*P_m.*(1-vol))./(vol+eta.*(1-vol));
K_T = 1./P;
m=1+4*K_T.*ni12.^2./e11;
e22 = (4*K_T.*g23)./(K_T+m.*g23);
e22(vol==0)=e11_m; e22(vol==1)=e11_f;
e33 = e22;
    
ni23 = ni23_f.*vol+ni23_m.*(1-vol).*(1+ni12_m-ni12.*(e11_m./e11))./...
     (1-ni12_m.^2+ni12_m.*ni12.*(e11_m./e11));
ni32 = e33./e22.*ni23; ni21 = e22./e11.*ni12;
ni31 = e33./e11.*ni13; 
g23 = e22./2./(1+ni23);

delta = 1-ni12.*ni21-ni23.*ni32-ni31.*ni13-2.*ni21.*ni32.*ni13;
q11 = e11.*(1-ni23.*ni32)./delta;
q12 = e11.*(ni21+ni31.*ni23)./delta;
q13 = e11.*(ni31+ni21.*ni32)./delta;
q22 = e22.*(1-ni31.*ni13)./delta;
q23 = e22.*(ni32+ni12.*ni31)./delta;
q33 = e33.*(1-ni12.*ni21)./delta;
q44 = g23; q55=g13; q66=g12;


%---------------------- END OF CODE---------------------- 

% ================ [functionalHahn.m] ================  
