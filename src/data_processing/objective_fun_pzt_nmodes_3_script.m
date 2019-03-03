% call objective_fun for test cases calculated in SASE models

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;

D0 = 10e3;%  cut off frequency for high-pass Butterworth filter, double , Units: [Hz]
Nb=1; % Butterworth filter order, integer
w=12; % window size in points
overwrite = false; % allow overwriting existing results if true
number_of_modes_considered = 3; % number of modes considered in calculation of objective function score

% create path to the experimental processed data folder
data_path = fullfile( projectroot, 'data','raw','exp', 'pzt', filesep );

% create path to the numerical model data folder
modelfolder = 'SASE';
modelname = {'SASE1'}; % cell list of models for processing
%modelname = {'SASE1','SASE2','SASE9','SASE10'}; % cell list of models for processing

% cell lists of corresponding experimental files
% signals measured in pzt sensors at angles [0:15:90] in respect to excitation
exp_filename = {'pzt_simul_289x289p_HANN100_x30_10Vpp_200Hz',...
                'pzt_simul_289x289p_HANN100_x30_10Vpp_200Hz',...
                'pzt_simul_289x289p_HANN100_x30_10Vpp_200Hz',...
                'pzt_simul_289x289p_HANN100_x30_10Vpp_200Hz'}; 
% filename of interim data
interim_filename = ['objective_fun_simul_pzt_score_nmodes_',num2str(number_of_modes_considered)];

for iModel = 1:length(modelname)
    
    % prepare model input and output path
    model_output_path = prepare_model_paths('interim','num',modelfolder,modelname{iModel});
    model_input_path = prepare_model_paths('raw','num',modelfolder,modelname{iModel});
    folder  = model_input_path;
    list    = dir(fullfile(folder, '*input.mat')); % list of mat files to be processed
    nFile   = length(list);
    success = false(1, nFile);
    objective_fun_pzt_score=zeros(nFile,1);
    % check if output already exist
    if(overwrite||(~overwrite && ~exist([model_output_path,filesep,interim_filename,'.mat'], 'file')))
        try
            disp('calculating objective function score ...');
            % load experimental data file (signals)
            load([data_path,exp_filename{iModel}]); % signals, beta, L, time 
            for test_case = 1:nFile
                 %% START DATA PROCESSING              
                output_name = [model_input_path,filesep,num2str(test_case),'output'];
                % load numerical data file
                load(output_name); % FREQ CG wavenumber
                [score] = objective_fun_pzt(time,signals(:,1:512),L,CG,FREQ,wavenumber,number_of_modes_considered,w,D0,Nb);
                objective_fun_pzt_score(test_case,1)=score;
                %% END DATA PROCESSING
            end
            fprintf('Model: %s successfully processed\n', modelname{iModel});
            % save processed data to interim (intermidiate) data folder
            save([model_output_path,filesep,interim_filename],'objective_fun_pzt_score');
        catch
            fprintf('Failed model: %s\n', modelname{iModel});
        end
    else
        fprintf('Filename: %s for model\n%s already exist\n', interim_filename,modelname{iModel});
    end
end