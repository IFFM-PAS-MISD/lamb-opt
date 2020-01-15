function [ObjV] = obj_ga_plain_weave_known_mass(Phen,Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,h_p,h_f,h_w,a_f,a_w,g_f,g_w,fiberType,rho,a,b)
% OBJ_GA_UNIDIRECTIONAL   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = obj_ga_plain_weave(input1,input2,input3) 
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
%run([src_path,filesep,'models',filesep,'SASE',filesep,'inputs',filesep,'Fabric_1.m']);
%run([src_path,filesep,'models',filesep,'SASE',filesep,'inputs',filesep,'Fabric_2.m']);
%run([src_path,filesep,'models',filesep,'SASE',filesep,'inputs',filesep,'Fabric_3.m']);
%run([src_path,filesep,'models',filesep,'SASE',filesep,'inputs',filesep,'Fabric_4.m']);
%run([src_path,filesep,'models',filesep,'SASE',filesep,'inputs',filesep,'Fabric_6.m']);
ObjV = zeros(size(Phen,1),1);
fiberType = repmat(fiberType,[size(Phen,1),1]);
h_p = repmat(h_p,[size(Phen,1),1]);
h_f = repmat(h_f,[size(Phen,1),1]);
h_w = repmat(h_w,[size(Phen,1),1]);
a_f = repmat(a_f,[size(Phen,1),1]);
a_w = repmat(a_w,[size(Phen,1),1]);
g_f = repmat(g_f,[size(Phen,1),1]);
g_w = repmat(g_w,[size(Phen,1),1]);

%for k=1:size(Phen,1)
parfor k=1:size(Phen,1)
    %[k,size(Phen,1)]
    rho_m = Phen(k,1);
    e11_m = Phen(k,2)/1e9;
    e11_f = Phen(k,3)/1e9;
    ni12_m = Phen(k,4);
    ni12_f = Phen(k,5);
    vol_0 = Phen(k,6);
    rho_f = (rho - rho_m*(1-vol_0))/vol_0;
    e22_f = 0.1*e11_f;
    ni23_f =  ni12_f ;
    %% Mechanical properties  
    
     [Q11,Q12,Q13,Q21,Q22,Q23,Q31,Q32,Q33,Q44,Q55,Q66,rho1] = ...
            compfabricprop(fiberType(k,:), h_p(k), h_f(k), h_w(k), a_f(k), a_w(k), g_f(k), g_w(k), vol_0, ...
            e11_m, ni12_m, rho_m, e11_f, e22_f, ni12_f, ni23_f, rho_f,false);
%     [Q11,Q12,Q13,Q21,Q22,Q23,Q31,Q32,Q33,Q44,Q55,Q66,rho] = ...
%         compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_0, ...
%         e11_m, ni12_m, rho_m, e11_f, e22_f, ni12_f, ni23_f, rho_f,false);
        
    %% SASE
    [wavenumber,CG,FREQ] = main_SASE(rho1,Q11,Q12,Q13,Q22,Q23,Q33,Q44,Q55,Q66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
    %[rho rho1]
    [score] = objective_fun_mod(Data_polar,fmax,FREQ,number_of_modes_considered);
    %ObjV(k)=1-score;
    %ObjV(k)=0.5-score;
    ObjV(k)=a*(-1)*score+b;
end

%---------------------- END OF CODE---------------------- 

% ================ [obj_ga_plain_weave.m] ================  
