clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;

% create path to the experimental raw data folder
interim_data_path=fullfile( projectroot, 'data','interim','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
filename = 'interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_'; 


% load interim experimental data file
disp('loading data ...');
load([interim_data_path,filename]); % KXKYF_ kx_vec ky_vec f_vec

Data=KXKYF_; clear KXKYF_;
beta = 0:15:90;
b = beta*pi/180;
number_of_angles = length(beta);
fmax=f_vec(end);
%% check NAN
[m1,n1,number_of_frames]=size(Data); % number_of_frames is equal to number_of_frequencies
for i=1:m1
    for j=1:n1
        for k=1:number_of_frames
            if(isnan(Data(i,j,k)))
                Data(i,j,k)=0;
            end
        end
    end
end

%%
% input
% lxmax=kx_vec(end); % length
% lymax=ky_vec(end); % width
lxmax=kx_vec(end)/2; % length
lymax=ky_vec(end)/2; % width
lxmin=kx_vec(floor(length(kx_vec)/2)+1); % lxmin=0; first quarter
lymin=ky_vec(floor(length(kx_vec)/2)+1); % lymin=0; first quarter
% Define the resolution of the grid:
N=max([m1,n1]); % # no of grid points for R coordinate
if(mod(N,2)) N=N-1; end;
%%
disp('Preliminary calculation...');
% Polar data allocation: angle, radius(wavenumbers), time(frequency)
Data_polar=zeros(number_of_angles,N,number_of_frames);
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
    wavenumber_step(k)=(wavenumber_max(k)-wavenumber_min(k))/(N-1); % wavenumber step [1/m]
end
x=zeros(number_of_angles,N);
y=zeros(number_of_angles,N);
for k=1:number_of_angles 
    R=linspace(wavenumber_min(k),wavenumber_max(k),N);
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
[m1,n1,p1]=size(Data_polar);
for i=1:m1
    for j=1:n1
        for k=1:p1
            if(isnan(Data_polar(i,j,k)))
                Data_polar(i,j,k)=0;
            end
        end
    end
end
           
%% save data in polar coordinate system
disp('Saving data...');
% create path to the experimental processed data folder
processed_data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of processed data
processed_filename = ['polar_',filename];

% save processed data to processed data folder
save([processed_data_path,processed_filename],'Data_polar','x','y','wavenumber_max','fmax','beta','-v7.3');
param_filename = [processed_filename,'_param'];
number_of_wavenumber_points = N;
save([processed_data_path,param_filename],'x','y','wavenumber_max','fmax','beta','number_of_wavenumber_points');

