% This script extracts signals from polar wavefield
% at angles beta and selected radius
% which simulates pzt sensors 
% the results are stored in experimental raw data
% pzt subfolder

clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;
% input
overwrite = false; % allow overwriting existing results if true
beta = 0:15:90;
R=0.14; % radius of pzt placement
test_case=[14]; % select file numbers for processing
foldername='pzt'; % subfolder for storing results
% create path to the experimental input data folder
input_data_path = fullfile( projectroot, 'data','processed','exp','wavefield', filesep);

% create path to the experimental raw data folder
processed_data_path = prepare_folder_paths('raw','exp',foldername);

% filenames of data to be processed
% full field measurements
list = {'101x101p_Chirp10_0-250kHz', ...                % 1                           
        '251x251p_0-1200kHz_Chirp_10us_x100', ...       % 2  
        '251x251p_400kHz_10HC_x10', ...                 % 3 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_200Hz', ...  % 4 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_5000Hz', ... % 5 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_500Hz', ...  % 6 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_50Hz', ...   % 7 
        '289x289p_0-350kHz_CHP160_x3_2Vpp_200Hz', ...   % 8
        '289x289p_0-350kHz_CHP20_x3_18Vpp_50Hz', ...    % 9
        '289x289p_0-350kHz_CHP40_x3_18Vpp_50Hz', ...    % 10 
        '289x289p_0-350kHz_CHP80_x3_18Vpp_50Hz', ...    % 11 
        '289x289p_HANN500_x3_18Vpp_50Hz', ...           % 12 
        '289x289p_10HC400kHz_x3_7Vpp_50Hz', ...         % 13
        '289x289p_HANN100_x30_10Vpp_200Hz', ...         % 14 
        '289x289p_HANN100_x3_10Vpp_200Hz', ...          % 15 
        '289x289p_HANN25_x3_10Vpp_200Hz', ...           % 16 
        '289x289p_HANN50_x3_10Vpp_200Hz', ...           % 17 
        '289x289p_HANN50_x3_15Vpp_200Hz', ...           % 18 
        '289x289p_HANN50_x3_20Vpp_200Hz', ...           % 19 
        '289x289p_HANN50_x3_5Vpp_200Hz', ...            % 20 
        '493x493p_0-350kHz_CHP160_x3_18Vpp_50Hz', ...   % 21 
        '493x493p_HANN100_x10_10Vpp_200Hz', ...         % 22 
        '493z493p_0-350kHz_CHP20N_x3', ...              % 23
        '493z493p_0-350kHz_CHP40N_x3', ...              % 24 
        '493z493p_0-500kHz_CHP_40N_x3'};                % 25

    
disp('Extracting signals from polar wavefield');

nFile   = length(test_case);
success = false(1, nFile);
for k = test_case
    filename = list{k};
    processed_filename = ['pzt_simul_',filename]; % filename of processed .mat data
    % check if already exist
    if(overwrite||(~overwrite && ~exist([processed_data_path,processed_filename,'.mat'], 'file')))
        try 
            % load interim experimental data file
            disp('loading data');
            load([input_data_path,'polar_wavefield_',filename]); % Data_polar, number_of_points, radius, time, beta
            %% PROCESS DATA
            fprintf('Processing:\n%s\n',filename);
            Rpoints=linspace(0,radius,number_of_points);
            [A,I] = min(abs(Rpoints-R));
            Actual_radius = Rpoints(I);
            signals = squeeze(Data_polar(:,I,:));
            L=Actual_radius;
            %% END OF PROCESSING
            % save data in polar coordinate system
            disp('Saving data...');
            % save processed data to output data folder
            save([processed_data_path,processed_filename],'signals','time','L','beta');
            fprintf('Successfully processed:\n%s\n', filename);% successfully processed
        catch
            fprintf('Failed: %s\n', filename);
        end
    else
        fprintf('Filename: \n%s \nalready exist\n', processed_filename);
    end
end
