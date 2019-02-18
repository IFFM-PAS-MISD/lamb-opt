% This script transforms wavenumber-frequency data from Cartesian coordinates
% to polar coordinates at selected angles beta by using interpolation
% all *.mat files in interim experimental data folder are processed 
% and the results are stored in experimental processed data folder

clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;
% input
overwrite = false; % allow overwriting existing results if true
beta = 0:15:90;

% create path to the experimental interim data folder
interim_data_path = fullfile( projectroot, 'data','interim','exp', filesep );

% filenames of data to be processed
% full field measurements after 3D FFT transform (positive quarter)

% create path to the experimental processed data folder
processed_data_path = fullfile( projectroot, 'data','processed','exp', filesep );
disp('Cartesian to polar transform ');
folder  = interim_data_path;
list    = dir(fullfile(folder, '*.mat')); % list of mat files to be processed
nFile   = length(list);
success = false(1, nFile);
for k = 1:nFile
    filename = list(k).name;
    processed_filename = ['polar_',filename]; % filename of processed .mat data
    % check if already exist
    if(overwrite||(~overwrite && ~exist([processed_data_path,processed_filename], 'file')))
        try 
            % load interim experimental data file
            disp('loading data');
            load([interim_data_path,filename]); % KXKYF_ kx_vec ky_vec f_vec

            % input
            fmax=f_vec(end);
            kxmax=kx_vec(end);
            kymax=ky_vec(end);
            % end of input

            %% PROCESS DATA
            fprintf('Processing:\n%s\n',filename);
            [Data_polar,number_of_wavenumber_points,wavenumber_max] = cartesian_to_polar_wavefield(KXKYF_,kxmax,kymax,beta); 

            %% END OF PROCESSING
            % save data in polar coordinate system
            disp('Saving data...');
            % save processed data to processed data folder
            save([processed_data_path,processed_filename],'Data_polar','wavenumber_max','fmax','beta','-v7.3');
            [filepath,name,ext] = fileparts(filename);
            param_filename = ['polar_',name,'_param'];
            save([processed_data_path,param_filename],'wavenumber_max','fmax','beta','number_of_wavenumber_points');
            fprintf('Successfully processed:\n%s\n', filename);% successfully processed
        catch
            fprintf('Failed: %s\n', filename);
        end
    else
        fprintf('Filename: \n%s \nalready exist\n', processed_filename);
    end
end
