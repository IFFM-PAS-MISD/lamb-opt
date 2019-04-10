function ownerElement = findOwnerElement(Nx,Ny,nodeCoordinates_sem,elementNodes_sem,shape_order,Eps)
% FINDOWNERELEMENT   One line description of what the function or script performs (H1 line) 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = findOwnerElement(input1,input2,input3) 
% 
% Inputs: 
%    Nx,Ny - Number of points in uniform mesh in x,y direction, integer 
%    nodeCoordinates_sem - Description, logical, dimensions [m, n], Units: m 
%    elementNodes_sem - Spectral nodes topology, double, dimensions [m, n], Units: N
%    shape_order - element approximation order, integer (N=3,4,5,6,7,8,9)
%    Eps - accuracy, double, dimensions [m, n], Units: -
% 
% Outputs: 
%    ownerElement - a vector of the numbers of the elements
%                   containing the points of the regular grid,
%                   integer, dimensions Nx*Ny, Units: - 
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
L = max(nodeCoordinates_sem(:,1))-min(nodeCoordinates_sem(:,1))-2*Eps;
xi = min(nodeCoordinates_sem(:,1))+Eps:L/(Nx-1):max(nodeCoordinates_sem(:,1))-Eps;
B = max(nodeCoordinates_sem(:,2))-min(nodeCoordinates_sem(:,2))-2*Eps;
yi = min(nodeCoordinates_sem(:,2))+Eps:B/(Ny-1):max(nodeCoordinates_sem(:,2))-Eps;
[XI,YI] = meshgrid(xi,yi);
xp = reshape(XI,1,Nx*Ny);
yp = reshape(YI,1,Nx*Ny);
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
