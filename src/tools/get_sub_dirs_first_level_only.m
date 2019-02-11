function [sub_dirs_names] = get_sub_dirs_first_level_only(parent_dir)
% GET_SUB_DIRS_FIRST_LEVEL_ONLY   Get a list of all files and folders in this folder 
% 
% Syntax: [sub_dirs_names] = get_sub_dirs_first_level_only(parent_dir) 
% 
% Inputs: 
%    parent_dir - parent directory full path, string 
% 
% Outputs: 
%    sub_dirs_names - list of subdirectory folders and files, cell of strings 
% 
% Example: 
%    [sub_dirs_names] = get_sub_dirs_first_level_only(parent_dir) 
%    [sub_dirs_names] = get_sub_dirs_first_level_only('E:\work\projects\opus15\lamb-opt\reports\figures\')  
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also:  
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

files    = dir(parent_dir);
names    = {files.name};
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir] & ~strcmp(names, '.') & ~strcmp(names, '..');
% Extract only those that are directories.
sub_dirs_names = names(dirFlags);

%---------------------- END OF CODE---------------------- 

% ================ [get_sub_dirs_first_level_only.m] ================  
