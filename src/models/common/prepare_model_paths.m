function model_output_path = prepare_model_paths(modelname,data_process_type,data_origin)
% PREPARE_MODEL_PATHS   function creates model_output_path 
%     Also set current directory in folder related to the model
% Syntax: model_output_path = prepare_model_paths(modelname,data_process_type,data_origin)
% 
% Inputs: 
%    modelname -  name of model, string  
%    data_process_type - string: 'raw', 'interim', 'processed' 
%    data_origin - string: 'exp', 'num' 
% 
% Outputs:
%    model_output_path - path to store the result, string
%
% Example: 
%    model_output_path = prepare_model_paths(modelname,'raw','num'); 
%    model_output_path = prepare_model_paths('Model2','raw','num'); 
% 
% Other m-files required: none 
% Subfunctions: CHECK_INPUTS_OUTPUTS 
% MAT-files required: none 
% See also: 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

% load projectroot path
load project_paths projectroot src_path;

% set current path to model folder
model_path = [src_path,'models',filesep,modelname,filesep];
cd(model_path);

% create path to the numerical raw data folder (absolute path)
raw_data_path = fullfile( projectroot, 'data',data_process_type,data_origin, filesep );
model_output_path = [raw_data_path,modelname,'_out'];

% check if folder exist, if not create it
if ~exist(model_output_path, 'dir')
    mkdir(model_output_path);
end
end

%---------------------- END OF CODE---------------------- 

% ================ [prepare_model_paths.m] ================  