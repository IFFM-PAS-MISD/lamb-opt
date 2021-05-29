function [ObjV] = obj_ga_pzt_unidirectional_A0(Phen,time,signals,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,number_of_modes_considered,rho,w,D0,Nb,L)
% OBJ_GA_C_TENSOR   One line description of what the function or script performs (H1 line) 
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
ObjV = zeros(size(Phen,1),1);

% fittnes function scaling factors
a=100;
b=250;
D1 = 500e3;%
%for k=1:size(Phen,1)
signal1=squeeze(signals(:,:,1));
signal2=squeeze(signals(:,:,2));
w1=w(1);
w2=w(2);
D0_1=D0(1);
D0_2=D0(2);
Nb1=Nb(1);
Nb2=Nb(2);
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
    [wavenumber,CG,FREQ] = main_SASE(rho,Q11,Q12,Q13,Q22,Q23,Q33,Q44,Q55,Q66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);

    [score1,~] = objective_fun_pzt2(time,signal1,L,CG,FREQ,wavenumber,number_of_modes_considered,w1,D0_1,Nb1);
    [score2,~] = objective_fun_pzt2(time,signal2,L,CG,FREQ,wavenumber,number_of_modes_considered,w2,D0_2,Nb2);
    %[score] = objective_fun_pzt_selected_mode2(time,signals,L,CG,FREQ,wavenumber,[1,1,1,1,1,1,1],w,D0,D1,Nb)
    score=score1+10*score2;
    ObjV(k)=a*(-1)*score+b;
end

%---------------------- END OF CODE---------------------- 

% ================ [obj_ga_.m] ================  
