function [Cnew] = transform_prop_3D(C,stack_dir,theta)
% TRANSFORM_PROP_3D   transforms material properties for given stacking direction and rotation angle 
% 
% Syntax: [output1,output2] = transform_prop_3D(input1,input2,input3) 
% 
% Inputs: 
%    C - matrix of elastic constants, double, dimensions [6, 6], Units: Pa 
%    stack_dir - stacking direction 1, 2, or 3, integer 
%    theta - rotation angle according to layup, double, Units: deg 
% 
% Outputs: 
%    Cnew - transformed matrix of elastic constants, double, dimensions [6, 6], Units: Pa 
% 
% Example: 
%    [Cnew] = transform_prop_3D(C,stack_dir,theta)
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: LAMINA_3D 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

theta = theta*pi/180;

if(abs(stack_dir)==1)
    tmparray = [1 2 3];
elseif(abs(stack_dir)==2)
    tmparray = [2 3 1];
else
    tmparray = [3 1 2];
end

if(stack_dir < 0)
    L(tmparray,1) = [0 -cos(theta) -sin(theta)]';
    L(tmparray,2) = [0 -sin(theta) +cos(theta)]';
    L(tmparray,3) = [-1 0 0]';
else
    L(tmparray,1) = [0 cos(theta) sin(theta)]';
    L(tmparray,2) = [0 -sin(theta) +cos(theta)]';
    L(tmparray,3) = [1 0 0]';
end

% transformation matrix
T(1,1) = L(1,1)^2;          T(1,2) = L(2,1)^2;          T(1,3) = L(3,1)^2;
T(1,4) = L(2,1)*L(3,1);     T(1,5) = L(1,1)*L(3,1);     T(1,6) = L(1,1)*L(2,1);

T(2,1) = L(1,2)^2;          T(2,2) = L(2,2)^2;          T(2,3) = L(3,2)^2;
T(2,4) = L(2,2)*L(3,2);     T(2,5) = L(1,2)*L(3,2);     T(2,6) = L(1,2)*L(2,2);

T(3,1) = L(1,3)^2;          T(3,2) = L(2,3)^2;          T(3,3) = L(3,3)^2;
T(3,4) = L(2,3)*L(3,3);     T(3,5) = L(1,3)*L(3,3);     T(3,6) = L(1,3)*L(2,3);

T(4,1) = 2*L(1,2)*L(1,3);   T(4,2) = 2*L(2,2)*L(2,3);   T(4,3) = 2*L(3,2)*L(3,3);
T(4,4) = L(2,2)*L(3,3)+L(2,3)*L(3,2);
T(4,5) = L(1,2)*L(3,3)+L(1,3)*L(3,2);
T(4,6) = L(1,2)*L(2,3)+L(1,3)*L(2,2);

T(5,1) = 2*L(1,1)*L(1,3);   T(5,2) = 2*L(2,1)*L(2,3);   T(5,3) = 2*L(3,1)*L(3,3);
T(5,4) = L(2,1)*L(3,3)+L(2,3)*L(3,1);
T(5,5) = L(1,1)*L(3,3)+L(1,3)*L(3,1);
T(5,6) = L(1,1)*L(2,3)+L(1,3)*L(2,1);

T(6,1) = 2*L(1,1)*L(1,2);   T(6,2) = 2*L(2,1)*L(2,2);   T(6,3) = 2*L(3,1)*L(3,2);
T(6,4) = L(2,1)*L(3,2)+L(2,2)*L(3,1);
T(6,5) = L(1,1)*L(3,2)+L(1,2)*L(3,1);
T(6,6) = L(1,1)*L(2,2)+L(1,2)*L(2,1);
      
Cnew = T'*C*T;

%---------------------- END OF CODE---------------------- 

% ================ [transform_prop_3D.m] ================  
