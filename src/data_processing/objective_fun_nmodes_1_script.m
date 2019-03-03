% call objective_fun for test cases calculated in SASE models

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;

overwrite = false; % allow overwriting existing results if true
number_of_modes_considered = 1; % number of modes considered in calculation of objective function score

% create path to the experimental processed data folder
data_path = fullfile( projectroot, 'data','processed','exp', filesep );

% create path to the numerical model data folder
modelfolder = 'SASE';
modelname = {'SASE1','SASE2','SASE9','SASE10'}; % cell list of models for processing
% cell lists of corresponding experimental files
% full field measurements after 3D FFT transform (1 quarter)
% and conversion to polar coordinate system
exp_filename = {'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF',...
                'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF',...
                'polar_interim_9x9p_HANN100_x30_10Vpp_200Hz_KXKYF_pzt_simul1',...
                'polar_interim_17x17p_HANN100_x30_10Vpp_200Hz_KXKYF_pzt_simul1'}; 
% filename of interim data
interim_filename = ['objective_fun_score_nmodes_',num2str(number_of_modes_considered)];

for iModel = 1:length(modelname)
    
    % prepare model input and output path
    model_output_path = prepare_model_paths('interim','num',modelfolder,modelname{iModel});
    model_input_path = prepare_model_paths('raw','num',modelfolder,modelname{iModel});
    folder  = model_input_path;
    list    = dir(fullfile(folder, '*input.mat')); % list of mat files to be processed
    nFile   = length(list);
    success = false(1, nFile);
    objective_fun_score=zeros(nFile,1);
    % check if output already exist
    if(overwrite||(~overwrite && ~exist([model_output_path,filesep,interim_filename,'.mat'], 'file')))
        try
            disp('calculating objective function score ...');
            % load experimental data file
            load([data_path,exp_filename{iModel}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points     
            for test_case = 1:nFile
                 %% START DATA PROCESSING
                output_name = [model_input_path,filesep,num2str(test_case),'output'];
                % load numerical data file
                load(output_name); % FREQ CG wavenumber
                [score] = objective_fun(Data_polar,fmax,FREQ,number_of_modes_considered);
                objective_fun_score(test_case,1)=score;         
                %% END DATA PROCESSING
            end
            fprintf('Model: %s successfully processed\n', modelname{iModel});
        catch
            fprintf('Failed model: %s\n', modelname{iModel});
        end
    else
        fprintf('Filename: %s for model\n%s already exist\n', interim_filename,modelname{iModel});
    end
        % save processed data to interim (intermidiate) data folder
        save([model_output_path,filesep,interim_filename],'objective_fun_score');
end