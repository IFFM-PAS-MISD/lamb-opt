function [obj_score] = objective_fun(Data_polar,fmax,FREQ)
% OBJECTIVE_FUN   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = objective_fun(input1,input2,input3) 
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
%    [output1,output2] = objective_fun(input1,input2,input3) 
%    [output1,output2] = objective_fun(input1,input2) 
%    [output1] = objective_fun(input1,input2,input3) 
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

[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);
obj_score=0;
for j=1:number_of_angles % beta
    %% Function to logical matrix
    % convert numerical dispersion curve into logical matrix

    H=logical(zeros(number_of_wavenumber_points,number_of_frequency_points));
    for i=1:number_of_wavenumber_points
        [~,I] = min(abs( FREQ(1,i,j) - fvec )); % mode 1
        H(i,I) = 1;
        [~,I] = min(abs( FREQ(2,i,j) - fvec )); % mode 2
        H(i,I) = 1;
        [~,I] = min(abs( FREQ(3,i,j) - fvec )); % mode 3
        H(i,I) = 1;
        [~,I] = min(abs( FREQ(4,i,j) - fvec )); % mode 4
        H(i,I) = 1;
    end
    start_idx1 = 4;
    start_idx2 = 2;
    end_idx1 = number_of_wavenumber_points -1;
    end_idx2 = number_of_frequency_points -1;
    % apply logical matrix to experimental dispersion curve matrix
    dispersion = H(start_idx1:end_idx1,start_idx2:end_idx2).*squeeze(abs(Data_polar(j,start_idx1:end_idx1,start_idx2:end_idx2)));
    
    [I,J] = find(H(start_idx1:end_idx1,start_idx2:end_idx2)==true);

    [m,n]=size(dispersion);
    % sum of values
    for k =1:length(J)
        obj_score = obj_score + ((dispersion(I(k),J(k))))/(m*n);
    end
end

%---------------------- END OF CODE---------------------- 

% ================ [objective_fun.m] ================  
