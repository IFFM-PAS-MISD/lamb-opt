function model1_test(test_case,model_output_path,overwrite)
% MODEL1-TEST computes elementwise square of matrix A defined in input 
% file no: test_case
% Inputs:
%       test_case: integer, case number corresponding to input number
%       model_output_path: string, path for saving the result
%       overwirite: logical, parameter defining if overwritting results is allowed  

% 
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
