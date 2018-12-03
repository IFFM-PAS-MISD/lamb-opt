function model1_test(test_case,model_output_path,overwrite)
% MODEL1_TEST   computes elementwise square of matrices 
%    matrices A, B are defined in input file  
%    filename is 'input'+integer definde by test_case 
% 
% Syntax: model1_test(test_case,model_output_path,overwrite)
% 
% Inputs: 
%    test_case - case number corresponding to input number, integer  
%    model_output_path - path for saving the result, string 
%    overwrite - parameter defining if overwritting results is allowed, logical 
% 
% Outputs: none
% 
% Example: 
%    model1_test(test_case,model_output_path,overwrite) 
%    model1_test(test_case,model_output_path,true) 
%    model1_test(10,'E:\data\raw\num\model1_out\outputs\output1',true) 
% 
% Other m-files required: none 
% Subfunctions: CHECK_INPUTS_OUTPUTS 
% MAT-files required: none 
% See also: 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

[input_filename,output_filename_marker,outputs_folder_case] = check_inputs_outputs(test_case,model_output_path);

% if overwrite == false check if this case is already calculated 
if ~overwrite
    % check if output file(s) for test_case exist
    if exist(output_filename_marker, 'file')
        disp(['Test case ',num2str(test_case),' already calculated']);
        return;
    end
end
disp(['Calculating test case ',num2str(test_case)]);
%% START FUNCTION CODE
% run script with input values
run(input_filename);
A = A.^2;
B = B.^2;
%% END FUNCTION CODE

% save values to output files
output1_filename = [outputs_folder_case,'out1.txt'];
output2_filename = [outputs_folder_case,'out2.txt'];
save(output1_filename,'A','-ascii');
save(output2_filename,'B','-ascii');
% add marker (empty file) indicating that the results for test_case exist
if ~exist(output_filename_marker, 'file' )
    fid = fopen(output_filename_marker,'w');
    fclose(fid);
end
disp(['Test case ',num2str(test_case),' finished']);
end

%---------------------- END OF CODE---------------------- 

% ================ [model1_test.m] ================  
