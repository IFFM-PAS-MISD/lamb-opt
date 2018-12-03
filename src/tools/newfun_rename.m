function newfun_rename(oldname,newname)
% NEWFUN_RENAME   is a function which renames a function 
%    Renaming includes headers and other commented locations
%    The old function is deleted
% 
% Syntax: newfun_rename(oldname,newname) 
% 
% Inputs: 
%    oldname - String, old function name 
%    newname - String, new function name 
% 
% Outputs: none
% 
% Example: 
%    newfun_rename('test_fun','renamed_test_fun') 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: NEWFUN 
%

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 
%
% Inspired by: Frank González-Morphy, frank.gonzalez-morphy@mathworks.de

%---------------------- BEGIN CODE---------------------- 

if nargin == 0, help(mfilename); return, end
if nargin == 1
    error('   MSG: You must enter at least the New Functionname!')
end

try
    newname = CheckIfFileExists(newname);
catch
    error('  MSG: NEWFCN_RENAME:CheckIfFileExists - After rechecking run again!')
end

ext_mfile = '.m';
oldname_ext = [oldname ext_mfile];
newname_ext = [newname ext_mfile];

% Open File for Read content
fid = fopen(oldname_ext, 'r');
if(fid==-1) 
    error(['  MSG: ',oldname_ext, '_file does not exist - After rechecking run again!'])
end
% Scaned_File = fread(fid, '*char')
Scaned_File   = fscanf(fid, '%c', inf);
fclose(fid);

% Upper Filenames
oldname_U = upper(oldname);
newname_U = upper(newname);

Modified_File = Scaned_File;
% replace strings
Modified_File = strrep(Modified_File,oldname,newname);
Modified_File = strrep(Modified_File,oldname_U,newname_U);

% Delete the old Function File
delete(oldname_ext);
% Open File for Writing 
fid = fopen(newname_ext, 'w+');
% Create and Write content into the File
fprintf(fid, '%s', Modified_File);
% fprintf(fid, '%c', Modified_File);
st = fclose(fid);
if st == 0
    disp(['   MSG: <' oldname '.m> renamed in :'])
    disp(['   MSG: <' newname '.m> successfuly renamed!'])
    edit(newname_ext);  % Edit the new Function on the Editor
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%    CheckIfFileExists    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function notexisting = CheckIfFileExists(fcnname)
% Sub-Function does check if File exist already, if so ask for overriding
% it or maybe save it under another name.

ex = exist(fcnname);  % does M-Function already exist ? Loop statement
while ex == 2         % rechecking existence
    overwrite = 0;    % Creation decision
    msg = sprintf(['Sorry, but Function -< %s.m >- does already exist!\n', ...
        'Do you wish to Overwrite it ?'], fcnname);
    % Action Question: Text, Title, Buttons and last one is the Default
    action = questdlg(msg, ' Overwrite Function?', 'Yes', 'No','No');
    if strcmp(action,'Yes') == 1
        ex = 0; % go out of While Loop, set breaking loop statement
    else
        % Dialog for new Functionname
        fcnname = char(inputdlg('Enter new Function Name ... ', 'NEWFCN - New Name'));
        if isempty(fcnname) == 1  % {} = Cancel Button => "1"
            error('   MSG: User decided to Cancel !')
        else
            ex = exist(fcnname);  % does new functionname exist ?
        end
    end
end
overwrite = 1;
if overwrite == 1
    notexisting = fcnname;
end

%---------------------- END OF CODE---------------------- 

% ================ [newfun_rename.m] ================  
