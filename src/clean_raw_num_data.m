% clean raw numerical data

% load projectroot path
load project_paths projectroot src_path;

disp('Removing raw numerical data');
% interim numerical data path
raw_num_path=fullfile( projectroot, 'data','raw','num',filesep );
% Get a list of all folders in this folder
[sub_dirs_names] = get_sub_dirs_first_level_only(raw_num_path);
% remove all folders and files in figures folder
if(isempty(sub_dirs_names))
    disp('Raw numerical data folder is already empty')
else
    for k = 1 : length(sub_dirs_names)
        fprintf('Removing folder #%d = %s\n', k, sub_dirs_names{k});
        rmdir([raw_num_path,sub_dirs_names{k}],'s');
    end
end
