% clean all data except raw data (subfolders only)

% load projectroot path
load project_paths projectroot src_path;

disp('Removing interim data');
% interim numerical data path
interim_num_path=fullfile( projectroot, 'data','interim','num',filesep );
% Get a list of all folders in this folder
[sub_dirs_names] = get_sub_dirs_first_level_only(interim_num_path);
% remove all folders and files in figures folder
if(isempty(sub_dirs_names))
    disp('Interim numerical data folder is already empty')
else
    for k = 1 : length(sub_dirs_names)
        fprintf('Removing folder #%d = %s\n', k, sub_dirs_names{k});
        rmdir([interim_num_path,sub_dirs_names{k}],'s');
    end
end

% interim experimental data path
interim_exp_path=fullfile( projectroot, 'data','interim','exp',filesep );
% Get a list of all folders in this folder
[sub_dirs_names] = get_sub_dirs_first_level_only(interim_exp_path);
% remove all folders and files in figures folder
if(isempty(sub_dirs_names))
    disp('Interim experimental data folder is already empty')
else
    for k = 1 : length(sub_dirs_names)
        fprintf('Removing folder #%d = %s\n', k, sub_dirs_names{k});
        rmdir([interim_exp_path,sub_dirs_names{k}],'s');
    end
end

disp('Removing processed data');
% processed numerical data path
processed_num_path=fullfile( projectroot, 'data','processed','num',filesep );
% Get a list of all folders in this folder
[sub_dirs_names] = get_sub_dirs_first_level_only(processed_num_path);
% remove all folders and files in figures folder
if(isempty(sub_dirs_names))
    disp('Processed numerical data folder is already empty')
else
    for k = 1 : length(sub_dirs_names)
        fprintf('Removing folder #%d = %s\n', k, sub_dirs_names{k});
        rmdir([processed_num_path,sub_dirs_names{k}],'s');
    end
end

% processed experimental data path
processed_exp_path=fullfile( projectroot, 'data','processed','exp',filesep );
% Get a list of all folders in this folder
[sub_dirs_names] = get_sub_dirs_first_level_only(processed_exp_path);
% remove all folders and files in figures folder
if(isempty(sub_dirs_names))
    disp('Processed experimental data folder is already empty')
else
    for k = 1 : length(sub_dirs_names)
        fprintf('Removing folder #%d = %s\n', k, sub_dirs_names{k});
        rmdir([processed_exp_path,sub_dirs_names{k}],'s');
    end
end
