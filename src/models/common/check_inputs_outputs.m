function [input_filename,output_filename_marker,outputs_folder_case] = check_inputs_outputs(test_case,model_output_path)
%CHECK_INPUTS_OUTPUTS - check if appropriate input and output
% directories and files exists, if not actions are taken
% Inputs:
%       test_case: integer, case number corresponding to input number
%       model_output_path: string, path for saving the result
% Outputs:
%       input_filename: string, relative path
%       output_filename_marker: string, '.exist' empty file
%       outputs_folder_case: string, output directory for test_case


% check if 'inputs' folder exist, if not create it
inputs_folder = ['.',filesep,'inputs',filesep];
if ~exist(inputs_folder, 'dir')
    mkdir(inputs_folder);
end
% check if 'outputs' folder exist, if not create it
outputs_folder = [model_output_path,filesep,'outputs',filesep];
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

