% This script transforms full wavefield data (2D-space-time)
% by using 3D FFT to wavenumber-frequency domain (kx-ky-f)
% only selected files (test_case)in raw experimental data folder are processed 
% and the results are stored in experimental interim data folder
% only positive quarter in wavenumber domain is stored

clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;
% input
overwrite = false; % allow overwriting existing results if true
% length and width of full wavefield area
Length = 0.726;    
Width = 0.726;
test_case=[1]; % select file numbers for processing

% create path to the experimental raw data folder
raw_data_path = fullfile( projectroot, 'data','raw','exp', filesep );

% create path to the experimental interim data folder
interim_data_path = fullfile( projectroot, 'data','interim','exp', filesep );

% filenames of data to be processed
% full field measurements
list = {'499x499p_chp200_x40_18Vpp_250Hz', ...          % 1                           
        '499x499p_chp200_x8_18Vpp_100Hz' };                % 2

disp('Spatial to wavenumber wavefield transform');
folder  = raw_data_path;
nFile   = length(test_case);
success = false(1, nFile);
for k = test_case
    filename = list{k};
    processed_filename = ['interim_',filename,'_KXKYF']; % filename of processed .mat data
    % check if already exist
    if(overwrite||(~overwrite && ~exist([interim_data_path,processed_filename,'.mat'], 'file')))
        try 
            % load raw experimental data file
            disp('loading data');
            load([raw_data_path,filename]); % Data Length Width time
            %% PROCESS DATA
            fprintf('Processing:\n%s\n',filename);
            [KXKYF_,kx_vec,ky_vec,f_vec] = spatial_to_wavenumber_wavefield(Data,Length,Width,time);
            %% END OF PROCESSING
            % save data in polar coordinate system
            disp('Saving data...');
            % save processed data to processed data folder
            save([interim_data_path,processed_filename],'KXKYF_','f_vec','kx_vec', 'ky_vec','-v7.3'); % use compression
            %save([interim_data_path,processed_filename],'KXKYF_','f_vec','kx_vec', 'ky_vec'); % no compression
            [filepath,name,ext] = fileparts(filename);
            fprintf('Successfully processed:\n%s\n', filename);% successfully processed
        catch
            fprintf('Failed: %s\n', filename);
        end
    else
        fprintf('Filename: \n%s \nalready exist\n', processed_filename);
    end
end
