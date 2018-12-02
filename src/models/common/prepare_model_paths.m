function model_output_path = prepare_model_paths(modelname,data_process_type,data_origin)
% PREPARE_MODEL_PATHS - function creates model_output_path
% Inputs:
%       modelname - string, name of model
%       data_process_type - string, 'raw', 'interim', 'processed' 
%       data_origin - string, 'exp', 'num'

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