function [ObjV] = obj_ga_C_tensor_known_mass_unidirectional_kx_ky(Phen,Data_polar,layup,h,fmin,fmax,wavenumber_max,number_of_frequency_points,beta,stack_dir,np,nele_layer,number_of_modes_considered,rho)
% obj_ga_C_tensor_known_mass_unidirectional_kx_ky   Objective function value in kx-ky plane
% 
% Syntax: [ObjV] = obj_ga_C_tensor_known_mass_unidirectional_kx_ky(Phen,Data_polar,layup,h,fmin,fmax,wavenumber_max,number_of_frequency_points,beta,stack_dir,np,nele_layer,number_of_modes_considered,rho)
% 
% Inputs: 
%    wavenumber - wavenumber matrix computed by SASE model for frequencies selected from experimental data, double [rad/m]
%    dimensions[number_of_modes,number_of_frequency_points,number_of_angles]
%    number_of_modes_considered - number of modes considered in calculation of the score
% 
% Outputs: 
%    output1 - Description, integer, dimensions [m, n], Units: - 
%    output2 - Description, double, dimensions [m, n], Units: m/s^2 
% 
% Example: 
%    [output1,output2] = obj_ga_plain_weave(input1,input2,input3) 
%    [output1,output2] = obj_ga_plain_weave(input1,input2) 
%    [output1] = obj_ga_plain_weave(input1,input2,input3) 
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

% load projectroot path
load project_paths projectroot src_path;
ObjV = zeros(size(Phen,1),1);

% fittnes function scaling factors
a=100;
b=250;
%for k=1:size(Phen,1)
parfor k=1:size(Phen,1)
    %[k,size(Phen,1)]
    Q11 = Phen(k,1);
    Q12 = Phen(k,2);
    Q13 = Phen(k,3);
    Q22 = Phen(k,4);
    Q23 = Phen(k,5);
    Q33 = Phen(k,6);
    Q44 = Phen(k,7);
    Q55 = Phen(k,8);
    Q66 = Phen(k,9);
      
    %% SASE
    
    [wavenumber,CG,FREQ] = main_SASE2(rho,Q11,Q12,Q13,Q22,Q23,Q33,Q44,Q55,Q66,layup,h,fmin,fmax,number_of_frequency_points,beta,stack_dir,np,nele_layer);
    [score] = objective_fun_kx_ky(Data_polar,wavenumber_max,wavenumber,number_of_modes_considered);
    %ObjV(k)=1-score;
    %ObjV(k)=0.5-score;
    ObjV(k)=a*(-1)*score+b;
end

%---------------------- END OF CODE---------------------- 

% ================ [obj_ga_plain_weave.m] ================  
