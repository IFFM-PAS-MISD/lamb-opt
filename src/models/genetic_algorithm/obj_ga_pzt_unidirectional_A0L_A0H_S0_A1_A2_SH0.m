function [ObjV] = obj_ga_pzt_unidirectional_A0L_A0H_S0_A1_A2_SH0(Phen,time,signals_A0L,signals_A0H,signals_S0,signals_A1,signals_A2,signals_SH0,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,number_of_modes,selected_mode_A0L,selected_mode_A0H,selected_mode_S0,selected_mode_A1,selected_mode_A2,selected_mode_SH0,rho,w_A0L,w_A0H,w_S0,w_A1,w_A2,w_SH0,D0_A0L,D0_A0H,D0_S0,D0_A1,D0_A2,D0_SH0,D1_A0L,D1_A0H,D1_S0,D1_A1,D1_A2,D1_SH0,Nb_A0L,Nb_A0H,Nb_S0,Nb_A1,Nb_A2,Nb_SH0,L)
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
    [wavenumber,~,FREQ] = main_SASE(rho,Q11,Q12,Q13,Q22,Q23,Q33,Q44,Q55,Q66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
    %% mode tracing automatic
    for j=1:7
        [FREQ_new] = mode_tracing_new(squeeze(FREQ(:,:,j)),wavenumber(:,j),number_of_modes);
        FREQ(1:number_of_modes,:,j)=FREQ_new';
    end
    FREQ=FREQ(1:number_of_modes,:,:);
    % A0 mode score at low frequency 
    [score_A0L] = objective_fun_pzt_selected_mode2(time,signals_A0L,L,FREQ,wavenumber,selected_mode_A0L,w_A0L,D0_A0L,D1_A0L,Nb_A0L);
    % A0 mode score at high frequency 
    [score_A0H] = objective_fun_pzt_selected_mode2(time,signals_A0H,L,FREQ,wavenumber,selected_mode_A0H,w_A0H,D0_A0H,D1_A0H,Nb_A0H);
    % S0 mode score
    [score_S0] = objective_fun_pzt_selected_mode2(time,signals_S0,L,FREQ,wavenumber,selected_mode_S0,w_S0,D0_S0,D1_S0,Nb_S0);
    % A1 mode score 
    [score_A1] = objective_fun_pzt_selected_mode2(time,signals_A1,L,FREQ,wavenumber,selected_mode_A1,w_A1,D0_A1,D1_A1,Nb_A1);
    % A2 mode score 
    [score_A2] = objective_fun_pzt_selected_mode2(time,signals_A2,L,FREQ,wavenumber,selected_mode_A2,w_A2,D0_A2,D1_A2,Nb_A2);
    % SH0 mode score 
    [score_SH0] = objective_fun_pzt_selected_mode2(time,signals_SH0,L,FREQ,wavenumber,selected_mode_SH0,w_SH0,D0_SH0,D1_SH0,Nb_SH0);
    
    score=score_A0L+50*score_A0H+10*score_S0+200*score_A1+100*score_A2+score_SH0;
    ObjV(k)=a*(-1)*score+b;
    [score_A0L,100*score_A0H,10*score_S0,200*score_A1,100*score_A2,score_SH0]
end

%---------------------- END OF CODE---------------------- 

% ================ [obj_ga_.m] ================  
