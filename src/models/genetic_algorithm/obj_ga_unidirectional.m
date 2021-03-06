function [ObjV] = obj_ga_unidirectional(Phen,Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered)
% OBJ_GA_UNIDIRECTIONAL   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = obj_ga_unidirectional(input1,input2,input3) 
% 
% Inputs: 
%    input1 - Description, string, dimensions [m, n], Units: ms 
%    input2 - Description, logical, dimensions [m, n], Units: m 
%    input3 - Description, double, dimensions [m, n], Units: N 
% 
% Outputs: 
%    output1 - Description, integer, dimensions [m, n], Units: - 
%    output2 - Description, double, dimensions [m, n], Units: m/s^2 
% 
% Example: 
%    [output1,output2] = obj_ga_unidirectional(input1,input2,input3) 
%    [output1,output2] = obj_ga_unidirectional(input1,input2) 
%    [output1] = obj_ga_unidirectional(input1,input2,input3) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

ObjV = zeros(size(Phen,1),1);
parfor k=1:size(Phen,1)
    %[k,size(Phen,1)]
    rhom = Phen(k,1);
    rhof = Phen(k,2);
    em = Phen(k,3);
    ef = Phen(k,4);
    nim = Phen(k,5);
    nif = Phen(k,6);
    vol = Phen(k,7);
    %% Mechanical properties  
    % homogenization of unidirectional fibre reinforce composite
    [rho,e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23] = homogenization(rhom,rhof,em,ef,nim,nif,vol);
    % Elastic constants of composite lamina in terms of principal material directions
    [C11,C12,C13,C22,C23,C33,C44,C55,C66] = lamina_3D(e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23);
    %% SASE
    [wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
    [score] = objective_fun(Data_polar,fmax,FREQ,number_of_modes_considered);
    %ObjV(k)=1-score;
    ObjV(k)=0.5-score;
end

%---------------------- END OF CODE---------------------- 

% ================ [obj_ga_unidirectional.m] ================  
