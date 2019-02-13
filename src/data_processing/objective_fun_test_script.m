% call objective_fun for test cases calculated in SASE1

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;

% create path to the numerical model data folder
foldername = 'SASE';
modelname = 'SASE1';

model_output_path = prepare_model_paths('raw','num',foldername,modelname);

% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
% and conversion to polar coordinate system
filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_'; 

% load experimental data file
disp('loading data ...');
load([data_path,filename]); % Data_polar x y wavenumber_max fmax beta number_of_wavenumber_points
number_of_modes_considered = 4; % number of modes considered in calculation of objective function score
%% START DATA PROCESSING
disp('calculating objective function score ...');
obj_score=zeros(121,1);
test_case = 0;
for i1 = 1:11
    for i2 = 1:11
        test_case= test_case+1;
        output_name = [model_output_path,filesep,num2str(test_case),'output'];
        % load numerical data file
        load(output_name); % FREQ CG wavenumber
        [score] = objective_fun(Data_polar,fmax,FREQ,number_of_modes_considered);
        obj_score(test_case,1)=score;
    end
end


%% END DATA PROCESSING

% create output path to the numerical interim data folder
interim_output_path = prepare_model_paths('interim','num',foldername,modelname);

% filename of processed data
interim_filename = 'objective_function_score';

% save processed data to interim (intermidiate) data folder
save([interim_output_path,filesep,interim_filename],'obj_score');

% add marker (empty file '.exist' ) indicating that 
% the processing results finished successfully and exist

output_filename_marker = [interim_output_path,filesep,'.exist'];

if ~exist(output_filename_marker, 'file' )
    fid = fopen(output_filename_marker,'w');
    fclose(fid);
end
