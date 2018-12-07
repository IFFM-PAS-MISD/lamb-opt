% convert function to logical matrix

clear all; close all;

model='SASE';
% load projectroot path
load project_paths projectroot src_path;

% create path to the numerical raw data folder
raw_data_path=fullfile( projectroot, 'data','raw','num',[model,'-out'], filesep );

% filename of data to be processed
filename = 'dispersion3D-3nodes-A0-4-2MHz.mat';

% load raw numerical data file
raw_data_num = load([raw_data_path,filename]);
%% START DATA PROCESSING

omega=raw_data_num.om(1:80:end);
wavenumber=raw_data_num.wavenumber(1:80:end);
plot(omega,wavenumber);
m=length(omega);
H=logical(sparse(m,m));
kw = max(wavenumber)*linspace(0,1,m); % liearly equally spaced vector of frequencies
%om = interp1(wavenumber,omega,kw','linear','extrap');
tic
for i=1:m
    [~,I] = min(abs(wavenumber(i)-kw));
    H(I,i) = 1;
end
toc
H = flipud(H); % only for checking figure; comment later on for algorithm
figure;
spy(H)

return;
%% END DATA PROCESSING

% create path to the experimental interim data folder
interim_data_path=fullfile( projectroot, 'data','interim','exp', filesep );

% filename of processed data
interim_filename = 'interim_data_testfile.txt';

% save processed data to interim (intermidiate) data folder
save([interim_data_path,interim_filename],'interim_data_exp','-ascii');

% add marker (empty file '.exist' ) indicating that 
% the processing results finished successfully and exist

output_filename_marker = [interim_data_path,'.exist'];

if ~exist(output_filename_marker, 'file' )
    fid = fopen(output_filename_marker,'w');
    fclose(fid);
end
