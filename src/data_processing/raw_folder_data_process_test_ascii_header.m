% process all files in raw data folder - test script template for ascii
% files without headers

clc;close all;clear all;

% allow overwriting existing results
overwrite=false;

% load projectroot path
load project_paths projectroot src_path;

% create path to the experimental raw data folder
raw_data_path=fullfile( projectroot, 'data','raw','exp', filesep );

% create path to the experimental interim data folder
interim_data_path=fullfile( projectroot, 'data','interim','exp', filesep );

folder  = raw_data_path;
list    = dir(fullfile(folder, '*.txt')); % list of mat files to be processed
nFile   = length(list);
success = false(1, nFile);
for k = 1:nFile
    file = list(k).name;
    interim_filename = ['interim_',file];
    exist_flag = exist( fullfile(interim_data_path,interim_filename), 'file'); % check if output file already exist
    
    if overwrite | (~overwrite & ~exist_flag) 
        try    
%% START DATA PROCESSING
%
        delimiterIn = ' '; % column separator
        headerlinesIn = 2; % number of header lines including colheaders
        A = importdata(fullfile(folder, file),delimiterIn,headerlinesIn);
        [m,n] = size(A.data);
        interim_data_exp = zeros(m,1);
        % compute RMS
        Data = A.data.^2; 
        for i=1:m
            interim_data_exp(i,1) = sqrt(sum(Data(i,:)))/n;
        end
        % save interim data
        save(fullfile(interim_data_path,interim_filename),'interim_data_exp','-ascii');
        fprintf('Success: %s\n', file);% succesfully processed
%% END DATA PROCESSING
        catch
            fprintf('Failed: %s\n', file);
        end          
    else
        disp(['File ',file,' has been already processed']);
    end  
end
