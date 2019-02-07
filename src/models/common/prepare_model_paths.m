function model_output_path = prepare_model_paths(data_process_type,data_origin,foldername,modelname)
% PREPARE_MODEL_PATHS   function creates model_output_path 
%     Also set current directory in folder related to the model
% Syntax: model_output_path = prepare_model_paths(modelname,data_process_type,data_origin)
% 
% Inputs:   
%    data_process_type - string: 'raw', 'interim', 'processed' 
%    data_origin - string: 'exp', 'num' 
%    foldername -  name of folder in which model resides, string
%    modelname -  name of model, string
% 
% Outputs:
%    model_output_path - path to store the result, string
%
% Example: 
%    model_output_path = prepare_model_paths('raw','num',modelname); 
%    model_output_path = prepare_model_paths('raw','num','Model2'); 
%    model_output_path = prepare_model_paths('raw','num','model','model1'); 
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

if nargin == 3
    modelname = foldername;
    % set current path to model folder
    model_path = [src_path,'models',filesep,modelname,filesep];
    cd(model_path);
end
if nargin == 4
    % set current path to model folder
    model_path = [src_path,'models',filesep,foldername,filesep];
    cd(model_path);
    modelname = [foldername,filesep,modelname];
end

% create path to the output data folder (absolute path)
data_output_path = fullfile( projectroot, 'data',data_process_type,data_origin, filesep );
model_output_path = [data_output_path,modelname,'_out'];

% check if folder exist, if not create it
if ~exist(model_output_path, 'dir')
    mkdir(model_output_path);
end
end

%---------------------- END OF CODE---------------------- 

% ================ [prepare_model_paths.m] ================  