clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;
% input
overwrite=false; % allow overwriting existing results if true
beta = 0:15:90;

% create path to the experimental raw data folder
interim_data_path=fullfile( projectroot, 'data','interim','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
filename = 'interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_'; 

% create path to the experimental processed data folder
processed_data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of processed data
processed_filename = ['polar_',filename];
if(overwrite||(~overwrite && ~exist([processed_data_path,processed_filename,'.mat'], 'file')))
    
    % load interim experimental data file
    disp('loading data');
    load([interim_data_path,filename]); % KXKYF_ kx_vec ky_vec f_vec

    % input
    Data=KXKYF_; clear KXKYF_;
    fmax=f_vec(end);
    kxmax=kx_vec(end)/2;
    kymax=ky_vec(end)/2;
    % end of input
    
    %% PROCESS DATA
    fprintf('Processing:\n %s\n', filename);
    [Data_polar,number_of_wavenumber_points,wavenumber_max] = cartesian_to_polar_wavefield(Data,kxmax,kymax,beta); 

    %% END OF PROCESSING
    % save data in polar coordinate system
    disp('Saving data...');

    % save processed data to processed data folder
    save([processed_data_path,processed_filename],'Data_polar','wavenumber_max','fmax','beta','-v7.3');
    param_filename = [processed_filename,'_param'];
    save([processed_data_path,param_filename],'wavenumber_max','fmax','beta','number_of_wavenumber_points');
else
    fprintf('Filename: \n%s \nalready exist\n', processed_filename);
end
