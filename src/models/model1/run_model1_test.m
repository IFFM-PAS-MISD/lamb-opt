% script to run model1 (template)

clc;close all;clear all;

% allow overwriting existing results
overwrite=false;

% retrieve model name based on foldername
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelname = pathstr(idx(end)+1:end);

% prepare model output path
model_output_path = prepare_model_paths(modelname,'raw','num');

%% START RUNNING MODEL
% replace 'model1_test' with your 'model_name'
for test_case = 1:2
    model1_test(test_case,model_output_path,overwrite);
end

%% END RUNNING MODEL
