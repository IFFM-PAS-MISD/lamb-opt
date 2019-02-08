% process raw data test script template

clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;

% create path to the experimental raw data folder
raw_data_path=fullfile( projectroot, 'data','raw','exp', filesep );

% filename of data to be processed
filename = '643p_pulse_5ms_x100_7Vpp_up250kHz.mat'; % measurements along a line


% load raw experimental data file
load([raw_data_path,filename]);

%% START DATA PROCESSING
%
Lenght = 0.9*sqrt(2);  

        % Unit vecotrs      
Fs =  1/(time(3)-time(2));    % sampling frequency
L = length(Data(1,:));
NFFT = 2^nextpow2(L); 
%f_vec = Fs/2*linspace(0,1,NFFT/2+1);  % frequency vector

[points frames] = size(Data);
deltay = Lenght/(points-1);                   
ky1 = (0:points-1)/points*2*pi/deltay;             % wavenumebr 1/m
%ky1 = (0:points-1)/points/deltay;             % angular wavenumebr rad/m
ky_vec = ky1-1/2*ky1(end);

% same as above
%deltaY = Lenght/(points-1);
%kY1 = mod( 1/2 + (0:(points-1))/points , 1 ) - 1/2;
%ky_vec2 = kY1*(2*pi/deltaY);


L_padded = 8*1024;
freqwave = fftshift(fftn(Data,[L_padded L_padded]));
f_vec = Fs/2*linspace(0,1,L_padded/2);  % frequency vector
%f_vec = Fs/2*linspace(0,1,L_padded/2)*2*pi;  % angular frequency vector
figure
set(gcf,'Color','w','Position',[100 100 640 480])
%imagesc(f_vec/10^3, ky_vec/1000, (abs(freqwave(:,floor(L_padded/2):L_padded))));   set(gca,'YDir','normal');  
C=(abs(freqwave(:,floor(L_padded/2):L_padded)));

%imagesc(f_vec, ky_vec, C);   set(gca,'YDir','normal');  %rad/m
[m,n]=size(C);
fmin=min(f_vec);
fmax=max(f_vec);
fvec2=fmin:(fmax-fmin)/(n-1):fmax;
kmin=min(ky_vec);
kmax=max(ky_vec);
kvec2=kmin:(kmax-kmin)/(m-1):kmax;
imagesc(fvec2, kvec2, C);   
set(gca,'YDir','normal');  
set(gca,'Fontsize',16)
%xlabel('Angular frequency (rad/s)')
xlabel('Frequency (Hz)')
ylabel('k_y (1/m)')
%colormap(cmap);
axis tight; 
b = max(caxis);  c = min(caxis);  
caxis([c b/5]); 
%interim_data_exp = 2*raw_data_exp;
return;
%% END DATA PROCESSING

% create path to the experimental interim data folder
interim_data_path=fullfile( projectroot, 'data','interim','exp', filesep );

% filename of processed data
interim_filename = 'interim_data_testfile.txt';

% save processed data to interim (intermidiate) data folder
save([interim_data_path,interim_filename],'interim_data_exp','-ascii');

% add marker (empty file '.exist' ) indicating that 
% the processing results finished successfully and exist

output_filename_marker = [interim_data_path,'.exist'];

if ~exist(output_filename_marker, 'file' )
    fid = fopen(output_filename_marker,'w');
    fclose(fid);
end
% comment
