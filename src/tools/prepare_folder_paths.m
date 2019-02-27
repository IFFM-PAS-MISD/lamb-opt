function output_path = prepare_folder_paths(data_process_type,data_origin,foldername)
% PREPARE_FOLDER_PATHS   function creates output_path 
%     Subfolder directory is created
% Syntax: output_path = prepare_folder_paths(modelname,data_process_type,data_origin)
% 
% Inputs:   
%    data_process_type - string: 'raw', 'interim', 'processed' 
%    data_origin - string: 'exp', 'num' 
%    foldername -  name of subfolder in which output will be stored, string
% 
% Outputs:
%    output_path - path to store the result, string
%
% Example: 
%    output_path = prepare_folder_paths('raw','num','foldername'); 
%    output_path = prepare_folder_paths('raw','num','wavefield'); 
% 
% Other m-files required: none 
% Subfunctions:  
% MAT-files required: none 
% See also: CHECK_INPUTS_OUTPUTS
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

% load projectroot path
load project_paths projectroot src_path;

% create path to the output data folder (absolute path)
data_output_path = fullfile( projectroot, 'data',data_process_type,data_origin, filesep );
output_path = [data_output_path,foldername,filesep];

% check if folder exist, if not create it
if ~exist(output_path, 'dir')
    mkdir(output_path);
end
end

%---------------------- END OF CODE---------------------- 

% ================ [prepare_folder_paths.m] ================  