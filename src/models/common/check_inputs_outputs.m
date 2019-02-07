function [input_filename,output_filename_marker,outputs_folder_case] = check_inputs_outputs(test_case,model_output_path)
% CHECK_INPUTS_OUTPUTS   Check exsitence of input and output folders/files
%    Check if appropriate input and output
%    directories and files exists, if not actions are taken 
% 
% Syntax: [input_filename,output_filename_marker,outputs_folder_case] = check_inputs_outputs(test_case,model_output_path)
% 
% Inputs: 
%    test_case - case number corresponding to input number, integer  
%    model_output_path - path for saving the result Description, string
% 
% Outputs: 
%    input_filename - relative path to input file, string 
%    output_filename_marker - path to '.exist' empty file including filename, string  
%    outputs_folder_case - output directory for test_case, string 
%
% Example: 
%    [input_filename,output_filename_marker,outputs_folder_case] = check_inputs_outputs(test_case,model_output_path) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also:  MODEL1-TEST
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

% check if 'inputs' folder exist, if not create it
inputs_folder = ['.',filesep,'inputs',filesep];
if ~exist(inputs_folder, 'dir')
    mkdir(inputs_folder);
end
% check if 'outputs' folder exist, if not create it
%outputs_folder = [model_output_path,filesep,'outputs',filesep];
outputs_folder = [model_output_path,filesep];
if ~exist(outputs_folder, 'dir')
    mkdir(outputs_folder);
end
% check if input file for test_case exist
input_filename = ['.',filesep,'inputs',filesep,'input',num2str(test_case),'.m'];
if ~exist(input_filename, 'file')
    error(['error: input',num2str(test_case),' does not exist']);
end
% check if output folder for test_case exist, if not create it
outputs_folder_case = [outputs_folder,'output',num2str(test_case),filesep];
if ~exist(outputs_folder_case, 'dir')
    mkdir(outputs_folder_case);
end
output_filename_marker = [outputs_folder_case,'.exist'];

end

%---------------------- END OF CODE---------------------- 

% ================ [check_inputs_outputs.m] ================  
