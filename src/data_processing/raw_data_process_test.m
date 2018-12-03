% process raw data test script template

clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;

% create path to the experimental raw data folder
raw_data_path=fullfile( projectroot, 'data','raw','exp', filesep );

% filename of data to be processed
filename = 'raw_data_testfile.txt';

% load raw experimental data file
raw_data_exp = load([raw_data_path,filename],'-ascii');

%% START DATA PROCESSING
%
interim_data_exp = 2*raw_data_exp;
%% END DATA PROCESSING

% create path to the experimental interim data folder
interim_data_path=fullfile( projectroot, 'data','interim','exp', filesep );

% filename of processed data
interim_filename = 'interim_data_testfile.txt';

% save processed data to interim (intermidiate) data folder
save([interim_data_path,interim_filename],'interim_data_exp','-ascii');

% add marker (empty file '.exist' ) indicating that the results exist
output_filename_marker = [interim_data_path,'.exist'];

if ~exist(output_filename_marker, 'file' )
    fid = fopen(output_filename_marker,'w');
    fclose(fid);
end
% comment
