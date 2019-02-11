% clean figures

% load projectroot path
load project_paths projectroot src_path;

% figure's path
figure_path=fullfile( projectroot, 'reports','figures', filesep );
% Get a list of all folders in this folder
[sub_dirs_names] = get_sub_dirs_first_level_only(figure_path);
% remove all folders and files in figures folder
disp('Removing figures folder');
if(isempty(sub_dirs_names))
    disp('Figures folder is already empty')
else
    for k = 1 : length(sub_dirs_names)
        fprintf('Removing folder #%d = %s\n', k, sub_dirs_names{k});
        rmdir([figure_path,sub_dirs_names{k}],'s');
    end
end


