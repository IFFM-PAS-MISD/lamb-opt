% script to run model1 (template)

clc;close all;clear all;

% allow overwriting existing results
overwrite=false;

% retrieve model name based on foldername
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelname = pathstr(idx(end)+1:end); % name of folder
%modelname = name; % name of running file

% prepare model output path
model_output_path = prepare_model_paths('raw','num',modelname);

%% START RUNNING MODEL
% replace 'model1_test' with your 'model_name'
for test_case = 1:2
    model1_test(test_case,model_output_path,overwrite);
end

%% END RUNNING MODEL
