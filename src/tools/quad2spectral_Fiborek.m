function [spec_element_nodes,spec_coords] = ...
    quad2spectral_Fiborek(elementNodes_fem,nodeCoordinates_fem,N)
% convert quad nodes mesh into spectral element mesh 
%    only for linear quad elements 
% 
% Syntax: [spec_element_nodes,spec_coords] = ...
%                                   quad2spec(elementNodes_fem,quad_coords,N) 
% 
% Inputs: 
%    elementNodes_fem - Quad nodes topology (element nodes), integer, dimensions [nQuadElements, 4]
%    nodeCoordinates_fem - coordinates of quad element nodes, double, dimensions [nQuadNodes, 3], Units: m 
%    N - element approximation order, integer (N=3,4,5,6,7,8,9)
% 
% Outputs: 
%    spec_element_nodes - spectral elements topology (element nodes), integer, dimensions [nQuadElements,(N+1)^2]
%    spec_coords - coordinates of spectral element nodes, integer, dimensions [nSpecNodes, 3], Units: m 
% 
%  13    14   15   16
%   O----O----O----O   
%   |    |    |    |
%  9O--10O----O11--O12
%   |    |    |    | 
%  5O---6O----O7---O8
%   |    |    |    |
%   O----O----O----O
%   1    2    3    4
%
% Example: 
%    [spec_element_nodes,spec_coords] = quad2spec(q,quad_coords,N)
% 
% Other m-files required: gll.m 
% Subfunctions: edgegeometry 
% MAT-files required: none 
% See also:
% 

% Author: Piotr Fiborek, M.Sc., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pfiborek@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 
n=N+1;
if n<3
    msgbox('You need at least 3 nodes on the edge of element', 'Error','error');
    return
elseif n>10
    msgbox('Maximum no of nodes on the edge of element is 10', 'Error','error');
    return
end
tic
numberElements = size(elementNodes_fem,1);
format long
n_int = n-2;
ksi = gll(n);
ksi_int = ksi(2:n-1);
% find mid point and length of every edge
[edgemidpointX edgelengthX] = edgegeometry(nodeCoordinates_fem,elementNodes_fem,1);
[edgemidpointY edgelengthY] = edgegeometry(nodeCoordinates_fem,elementNodes_fem,2);
[edgemidpointZ edgelengthZ] = edgegeometry(nodeCoordinates_fem,elementNodes_fem,3);

% the edge nodes coordinates
nodeCoordinatesX_temp = repmat(edgemidpointX,1,n_int)+...
    repmat(edgelengthX,1,n_int).*repmat(ksi_int,size(edgemidpointX,1),1);
nodeCoordinatesX_temp = reshape(nodeCoordinatesX_temp',[],1);
nodeCoordinatesY_temp = repmat(edgemidpointY,1,n_int)+...
    repmat(edgelengthY,1,n_int).*repmat(ksi_int,size(edgemidpointY,1),1);
nodeCoordinatesY_temp = reshape(nodeCoordinatesY_temp',[],1);
nodeCoordinatesZ_temp = repmat(edgemidpointZ,1,n_int)+...
    repmat(edgelengthZ,1,n_int).*repmat(ksi_int,size(edgemidpointZ,1),1);
nodeCoordinatesZ_temp = reshape(nodeCoordinatesZ_temp',[],1);

nodeCoordinates_int_edge = [nodeCoordinatesX_temp nodeCoordinatesY_temp nodeCoordinatesZ_temp];
[b1,~,n1] = unique(nodeCoordinates_int_edge,'first','rows');

nodeCoordinates_int_edge = b1;
unqNodes = 1:length(b1);
elementNodes_int_edge = reshape(unqNodes(n1),4*n_int,[])';
 
% coordinates of the nodes inside the element
elementNodes_int_int = (1:(numberElements*n_int.^2))+...
    max(max(elementNodes_int_edge));
elementNodes_int_int = reshape(elementNodes_int_int,n_int.^2,...
     numberElements);
elementNodes_int_int = elementNodes_int_int';
n_el1 = elementNodes_int_edge(:,1+0*n_int:1*n_int);
n_el2 = elementNodes_int_edge(:,1+1*n_int:2*n_int);
n_el3 = elementNodes_int_edge(:,3*n_int:-1:1+2*n_int);
n_el4 = elementNodes_int_edge(:,4*n_int:-1:1+3*n_int);
    
iSparse = repmat(1:numberElements,n_int,1);
jSparse = repmat((1:n_int:n_int*numberElements),n_int ,1)+repmat((0:n_int-1)',1,numberElements);
       
A13 = reshape(nodeCoordinates_int_edge(n_el3,2)-nodeCoordinates_int_edge(n_el1,2),[],n_int);
A13 = sparse(iSparse,jSparse,A13');
    
B13 = reshape(nodeCoordinates_int_edge(n_el1,1)-nodeCoordinates_int_edge(n_el3,1),[],n_int);
B13 = sparse(iSparse,jSparse,B13');
C13 = reshape((nodeCoordinates_int_edge(n_el1,2)-nodeCoordinates_int_edge(n_el3,2)).*...
        nodeCoordinates_int_edge(n_el1,1)+...
        (nodeCoordinates_int_edge(n_el3,1)-nodeCoordinates_int_edge(n_el1,1)).*...
        nodeCoordinates_int_edge(n_el1,2),[],n_int);
C13 = sparse(iSparse,jSparse,C13');
A42 = reshape(nodeCoordinates_int_edge(n_el2,2)-nodeCoordinates_int_edge(n_el4,2),[],n_int);
A42 = sparse(iSparse,jSparse,A42');
B42 = reshape(nodeCoordinates_int_edge(n_el4,1)-nodeCoordinates_int_edge(n_el2,1),[],n_int);
B42 = sparse(iSparse,jSparse,B42');
C42 = reshape((nodeCoordinates_int_edge(n_el4,2)-nodeCoordinates_int_edge(n_el2,2)).*...
        nodeCoordinates_int_edge(n_el4,1)+...
        (nodeCoordinates_int_edge(n_el2,1)-nodeCoordinates_int_edge(n_el4,1)).*...
        nodeCoordinates_int_edge(n_el4,2),[],n_int);
C42 = sparse(iSparse,jSparse,C42');
    
W = nonzeros(A13'*B42-B13'*A42);
Wx = nonzeros(B13'*C42-C13'*B42);
Wy = nonzeros(C13'*A42-A13'*C42);
    
X_int=reshape(Wx./W,[],1);
Y_int=reshape(Wy./W,[],1);
Z_int=ones(length(X_int),1).*unique(nodeCoordinates_int_edge(:,3));
nodeCoordinates_int_int = [X_int,Y_int,Z_int];

% regular spectral nodes topology
elementNodes_int = [elementNodes_int_edge,elementNodes_int_int]+double(max(max(elementNodes_fem)));
nodeCoordinates_int = [nodeCoordinates_int_edge;nodeCoordinates_int_int];
elementNodes_sem = [elementNodes_fem,elementNodes_int];
spec_coords = [nodeCoordinates_fem;nodeCoordinates_int];

elementNodes_pl3 = zeros(size(elementNodes_sem));
for ii = 1:n
    if ii == 1
        elementNodes_pl3(:,1+n*(ii-1):n*ii) = elementNodes_sem(:,[1,4+1:4+n_int,2]);
    elseif ii > 1 && ii < n
        elementNodes_pl3(:,1+n*(ii-1):n*ii) = elementNodes_sem(:,[4+4*n_int-ii+2,...
            (ii-2)*n_int+1+4+4*n_int:(ii-2)*n_int+4+4*n_int+n_int,...
            ii-1+4+1*n_int]);
    elseif ii == n
        elementNodes_pl3(:,1+n*(ii-1):n*ii) = elementNodes_sem(:,[4,...
            4+3*n_int:-1:4+2*n_int+1,3]);
    end
end
spec_element_nodes = elementNodes_pl3;

function [edgemidpoint, edgelength] = edgegeometry(nodeCoordinates,elementNodes,dim)
        
        nC = reshape(nodeCoordinates(elementNodes(:,[1,2,3 4 1]),dim),[],5);
        % midpoint of 4 edges for all elements
        edgemidpoint = (nC(:,2:5)+nC(:,1:4))/2;
        edgemidpoint = round(reshape(edgemidpoint',[],1)*1e6)*1e-6;
        
        % length of 4 edges for all elements
        edgelength = (nC(:,2:5)-nC(:,1:4))/2;
        edgelength = reshape(edgelength',[],1);

%---------------------- END OF CODE---------------------- 

% ================ [quad2spectral_Fiborek.m] ================  
