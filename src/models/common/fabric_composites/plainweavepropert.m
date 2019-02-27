function [Q11, Q12, Q13, Q21, Q22, Q23, Q31, Q32, Q33, Q44, Q55, Q66] = ...
    plainweavepropert(type, angle, q11, q12, q13, q22, q23, q33, q44, q55, q66)
% PLAINWEAVE_FILL transformation to global coord of effective properties of matrix, fill or warp) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = plainweavepropert(input1,input2,input3) 
% 
% Inputs: 
%    type - 'matrix', 'fill' or 'warp', string, Units:  
%    angle - fill or warp undulation angle for matrix zeros, dimensions [m, n], Units:
%    rad 
%    q11,...q66 - homogenized mechanical properties of matrix, fill or warp,
%               integer, dimensions [1, 1], Units: GPa
% 
% Outputs: 
%    Q11,...Q66 - global mechanical properties of matrix, 
%               fill or warp, integer, dimensions [m, n], Units: GPa 
% 
% Example: 
%    [output1,output2] = plainweavepropert(input1,input2,input3) 
%    [output1,output2] = plainweavepropert(input1,input2) 
%    [output1] = plainweavepropert(input1,input2,input3) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Piotr Fiborek, D.Sc., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pfiborek@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

switch type
    case 'fill'
        m = cos(angle);
        n = sin(angle);
        Q11 = q11.*m.^4+(2*q13+4*q55).*n.^2.*m.^2+q33.*n.^4;
        Q12 = q12.*m.^2+q23.*n.^2;
        Q13 = q13.*m.^4+(q11+(q33-4*q55)).*n.^2.*m.^2+q13.*n.^4;
        Q21 = q12.*m.^2+q23.*n.^2;
        Q22 = q22.*ones(size(angle));
        Q23 = q23.*m.^2+q12.*n.^2;
        Q31 = q13.*m.^4+(q11+(q33-4*q55)).*n.^2.*m.^2+q13.*n.^4;
        Q32 = q23.*m.^2+q12.*n.^2;
        Q33 = q33.*m.^4+(2*q13+4*q55).*n.^2.*m.^2+q11.*n.^4;
        Q44 = q44.*m.^2+q66.*n.^2;
        Q55 = q55.*m.^4+(q11+(-2*q13+(q33-2*q55))).*n.^2.*m.^2+q55.*n.^4;
        Q66 = q66.*m.^2+q44.*n.^2;
    case 'warp'
        m = cos(angle);
        n = sin(angle);
        Q11 = q22.*ones(size(angle));
        Q12 = q12.*m.^2+q23.*n.^2;
        Q13 = q23.*m.^2+q12.*n.^2;
        Q21 = q12.*m.^2+q23.*n.^2;
        Q22 = q11.*m.^4+(2*q13+4*q55).*n.^2.*m.^2+q33*n.^4;
        Q23 = q13.*m.^4+(q11+(q33-4*q55)).*n.^2.*m.^2+q13.*n.^4;
        Q31 = q23.*m.^2+q12.*n.^2;
        Q32 = q13.*m.^4+(q11+(q33-4*q55)).*n.^2.*m.^2+q13.*n.^4;
        Q33 = q33.*m.^4+(2*q13+4*q55).*n.^2.*m.^2+q11.*n.^4;
        Q44 = q55.*m.^4+(q11+(-2*q13+(q33-2*q55))).*n.^2.*m.^2+q55.*n.^4;
        Q55 = q44.*m.^2+q66.*n.^2;
        Q66 = q66.*m.^2+q44.*n.^2; 
    case 'matrix'
        Q11 = q11.*ones(size(angle));
        Q12 = q12.*ones(size(angle));
        Q13 = q13.*ones(size(angle));
        Q21 = q12.*ones(size(angle));
        Q22 = q22.*ones(size(angle));
        Q23 = q23.*ones(size(angle));
        Q31 = q13.*ones(size(angle));
        Q32 = q23.*ones(size(angle));
        Q33 = q33.*ones(size(angle));
        Q44 = q44.*ones(size(angle));
        Q55 = q55.*ones(size(angle));
        Q66 = q66.*ones(size(angle));
end

%---------------------- END OF CODE---------------------- 

% ================ [plainweavepropert.m] ================  
