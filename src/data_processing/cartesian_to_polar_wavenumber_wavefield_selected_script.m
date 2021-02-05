% This script transforms wavenumber-frequency data from Cartesian coordinates
% to polar coordinates at selected angles beta by using interpolation
% all *.mat files in interim experimental data folder are processed 
% and the results are stored in experimental processed data folder

clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;
% input
overwrite = false; % allow overwriting existing results if true
beta = 0:90/511:90;
selected_frequency_index = [40:40:300]; 
test_case=[1,3];
% create path to the experimental interim data folder
interim_data_path = fullfile( projectroot, 'data','interim','exp', filesep );

% filenames of data to be processed
% full field measurements after 3D FFT transform (positive quarter)
list = {'interim_499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni_KXKYF', ...          % 1  Length = 0.455;Width = 0.455;           
        'interim_499x499p_chp200_x40_6Vpp_250Hz_uni_KXKYF',... % 2 Length = 0.724;Width = 0.727;    
        'interim_483x483p_CHIRP_20-500kHz_125us_6Vpp_x3_stringer_intact_KXKYF'};                
% create path to the experimental processed data folder
processed_data_path = fullfile( projectroot, 'data','processed','exp', filesep );
disp('Cartesian to polar transform ');
folder  = interim_data_path;
nFile   = length(test_case);
for k = test_case
    filename = list{k};
    processed_filename = ['polar_f_',filename]; % filename of processed .mat data
    % check if already exist
    if(overwrite||(~overwrite && ~exist([processed_data_path,processed_filename,'.mat'], 'file')))
        try 
            % load interim experimental data file
            disp('loading data');
            load([interim_data_path,filename]); % KXKYF_ kx_vec ky_vec f_vec

            % input
            fmax=f_vec(end);
            kxmax=kx_vec(end);
            kymax=ky_vec(end);
            % end of input
            selected_frequencies = f_vec(selected_frequency_index);
            %% PROCESS DATA
            fprintf('Processing:\n%s\n',filename);
            [Data_polar,number_of_wavenumber_points,wavenumber_max] = cartesian_to_polar_wavenumber_wavefield_const_k(KXKYF_(:,:,selected_frequency_index),kxmax,kymax,beta); 

            %% END OF PROCESSING
            % save data in polar coordinate system
            disp('Saving data...');
            % save processed data to processed data folder
            save([processed_data_path,processed_filename],'Data_polar','wavenumber_max','wavenumber_max','selected_frequency_index','selected_frequencies','beta','-v7.3');
            [filepath,name,ext] = fileparts(filename);
            param_filename = ['polar_f_',name,'_param'];
            save([processed_data_path,param_filename],'wavenumber_max','wavenumber_max','selected_frequency_index','selected_frequencies','beta','number_of_wavenumber_points');
            fprintf('Successfully processed:\n%s\n', filename);% successfully processed
        catch
            fprintf('Failed: %s\n', filename);
        end
    else
        fprintf('Filename: \n%s \nalready exist\n', processed_filename);
    end
end
