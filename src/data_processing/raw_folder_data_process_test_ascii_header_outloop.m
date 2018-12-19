% process all files in raw data folder - test script template for ascii
% files without headers

clc;close all;clear all;

% allow overwriting existing results
overwrite=true;

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
Data = cell(nFile,1);
file_num = zeros(nFile,1);
cc = 0;
for k = 1:nFile
    file = list(k).name;
    interim_filename = ['interim_',file];
    exist_flag = exist( fullfile(interim_data_path,interim_filename), 'file'); % check if output file already exist
    
    if overwrite | (~overwrite & ~exist_flag) 
        try    
%% START LOADING DATA 
%
        delimiterIn = ' '; % column separator
        headerlinesIn = 2; % number of header lines including colheaders
        A = importdata(fullfile(folder, file),delimiterIn,headerlinesIn);
        temp =  A.data;
        success(k) = 1; % successfully loaded
%% END LOADING DATA
        catch
            fprintf('Failed: %s\n', file);
        end          
    else
        disp(['File ',file,' has been already processed']);
    end 
    if(success(k))
        cc=cc+1;
        Data{cc} = A.data; 
        file_num(cc) = k;
    end
end

%% START DATA PROCESSING
%   
ref_data = Data{1};
[m,n] = size(ref_data);
for k=2:cc
    interim_data_exp = Data{k} - ref_data; % make signal difference
    % save interim data
    file = list(file_num(k)).name;
    interim_filename = ['interim_',file];
    save(fullfile(interim_data_path,interim_filename),'interim_data_exp','-ascii');
    fprintf('Success: %s\n', file);% succesfully processed
    
end
%% END DATA PROCESSING
