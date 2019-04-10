function ownerElement = findOwnerElement(xp,yp,nodeCoordinates_sem,elementNodes_sem,shape_order)
% FINDOWNERELEMENT   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = findOwnerElement(input1,input2,input3) 
% 
% Inputs: 
%    xp,yp - points coordinates in uniform mesh in x,y direction, integer 
%    nodeCoordinates_sem - Description, logical, dimensions [m, n], Units: m 
%    elementNodes_sem - Spectral nodes topology, double, dimensions [m, n], Units: N
%    shape_order - element approximation order, integer (N=3,4,5,6,7,8,9)
%    % 
% Outputs: 
%    ownerElement - a vector of the numbers of the elements
%                   containing the points of the regular grid,
%                   integer, dimensions [length(xp)*length(yp),1], Units: - 
% 
% Example: 
%    [output1,output2] = findOwnerElement(input1,input2,input3) 
%    [output1,output2] = findOwnerElement(input1,input2) 
%    [output1] = findOwnerElement(input1,input2,input3) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: 
% 

% Author: Piotr Fiborek, M.Sc., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pfiborek@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 
disp('Determination of the owner element');
n = shape_order + 1;
elementNodes = elementNodes_sem(:,[1,n,n^2,n*(n-1)+1,1]);
ownerElement = zeros(length(xp),1);
for ipolygon = 1:size(elementNodes_sem,1)
    elementVertex = elementNodes(ipolygon,:);
    xv = nodeCoordinates_sem(elementVertex,1);
    yv = nodeCoordinates_sem(elementVertex,2);
    in = inpolygon(xp,yp,xv,yv);
    ownerElement(in) = ipolygon;
end

%---------------------- END OF CODE---------------------- 

% ================ [findOwnerElement.m] ================  
