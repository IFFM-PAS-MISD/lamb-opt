function [Data_polar,number_of_wavenumber_points,wavenumber_max] = cartesian_to_polar_wavefield(Data,kxmax,kymax,beta)
% CARTESIAN_TO_POLAR_WAVEFIELD   transform wavefield to polar coordinates 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
%    optional: more details about the function than in the H1 line 
% 
% Syntax: [output1,output2] = cartesian_to_polar_wavefield(input1,input2,input3) 
% 
% Inputs: 
%    input1 - Description, string, dimensions [m, n], Units: ms 
%    input2 - Description, logical, dimensions [m, n], Units: m 
%    input3 - Description, double, dimensions [m, n], Units: N 
% 
% Outputs: 
%    output1 - Description, integer, dimensions [m, n], Units: - 
%    output2 - Description, double, dimensions [m, n], Units: m/s^2 
% 
% Example: 
%    [output1,output2] = cartesian_to_polar_wavefield(input1,input2,input3) 
%    [output1,output2] = cartesian_to_polar_wavefield(input1,input2) 
%    [output1] = cartesian_to_polar_wavefield(input1,input2,input3) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

b = beta*pi/180;
number_of_angles = length(beta);

%% check NAN
[m1,n1,number_of_frames]=size(Data); % number_of_frames is equal to number_of_frequencies
% for i=1:m1
%     for j=1:n1
%         for k=1:number_of_frames
%             if(isnan(Data(i,j,k)))
%                 Data(i,j,k)=0;
%             end
%         end
%     end
% end
Data(isnan(Data))=0;
%%
% input
lxmax=kxmax; % length
lymax=kymax; % width
lxmin=0; % quarter
lymin=0; % quarter
% Define the resolution of the grid:
number_of_wavenumber_points=max([m1,n1]); % # no of grid points for R coordinate
if(mod(number_of_wavenumber_points,2)) % only even numbers
    number_of_wavenumber_points=number_of_wavenumber_points-1; 
end
%%
disp('Preliminary calculation...');
% Polar data allocation: angle, radius(wavenumbers), time(frequency)
Data_polar=zeros(number_of_angles,number_of_wavenumber_points,number_of_frames);
 %%
[XI,YI] = meshgrid(linspace(lxmin,lxmax,n1),linspace(lymin,lymax,m1)); % due to columnwise plotting n1 is for x coordinates and m1 is for y coordinates
X=reshape(XI,[],1);
Y=reshape(YI,[],1);
 %% 
wavenumber_max=zeros(number_of_angles ,1);
for k=1:number_of_angles 
    if(b(k) <= atan(lymax/lxmax))
        wavenumber_max(k,1)=sqrt((lxmax*tan(b(k))).^2+lxmax^2);
    else
        wavenumber_max(k,1)=sqrt((lymax*tan(pi/2-b(k))).^2+lymax^2);
    end
end
wavenumber_min = zeros(number_of_angles,1); % minimal wavenumbers [1/m]
wavenumber_step=zeros(number_of_angles,1);
for k=1:number_of_angles
    wavenumber_step(k)=(wavenumber_max(k)-wavenumber_min(k))/(number_of_wavenumber_points-1); % wavenumber step [1/m]
end
x=zeros(number_of_angles,number_of_wavenumber_points);
y=zeros(number_of_angles,number_of_wavenumber_points);
for k=1:number_of_angles 
    R=linspace(wavenumber_min(k),wavenumber_max(k),number_of_wavenumber_points);
    x(k,:) = R*cos(b(k));
    y(k,:) = R*sin(b(k));
end

 %%
 % convert Data from Cartesian to polar coordinates
 %%

 disp('Data conversion...');
% loop through time (frequency) frames
for frame=1:number_of_frames
    [frame,number_of_frames]
    ZI=Data(:,:,frame);
    Z=reshape(ZI,[],1);
    F = TriScatteredInterp(X,Y,Z,'linear');
    %Evaluate the interpolant at the locations (x, y).
    %The corresponding value at these locations is Ztemp:
    Zpolar = F(x,y);
    % store data
    Data_polar(:,:,frame)=Zpolar;
end

%% check NAN
%[m1,n1,p1]=size(Data_polar);
% for i=1:m1
%     for j=1:n1
%         for k=1:p1
%             if(isnan(Data_polar(i,j,k)))
%                 Data_polar(i,j,k)=0;
%             end
%         end
%     end
% end 
Data_polar(isnan(Data_polar))=0;
%---------------------- END OF CODE---------------------- 

% ================ [cartesian_to_polar_wavefield.m] ================  
