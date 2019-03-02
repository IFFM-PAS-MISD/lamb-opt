% dispersion compensation test
clc;close all;clear all;

% load projectroot path
load project_paths projectroot src_path;

% input
% central point
I=145; J=141;
% analyse Data(I,J+K,:) - K points on x axis from excitation
K=70;
dx=0.493/288;
L=K*dx; % dispersion compensation distance
D0 = 10e3;%  cut off frequency for high-pass Butterworth filter, double , Units: [Hz]
Nb=1; % Butterworth filter order, integer
w=12; % window size in points
overwrite = false; % allow overwriting existing results if true
% length and width of full wavefield area
Length = 0.493;    
Width = 0.493;
exp_test_case=14; % select file number from experiment for processing
num_test_case=121;%[25]; % select file numbers from SASE model for processing
modelfolder = 'SASE';
modelname = 'SASE1';
radians = false;
% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname,'_out'], filesep );

% create path to the experimental raw data folder
raw_data_path = fullfile( projectroot, 'data','raw','exp', filesep );

% create output path
output_path = prepare_model_paths('interim','exp',modelfolder,modelname);

% filenames of data to be processed
% full field measurements
list = {'101x101p_Chirp10_0-250kHz', ...                % 1                           
        '251x251p_0-1200kHz_Chirp_10us_x100', ...       % 2  
        '251x251p_400kHz_10HC_x10', ...                 % 3 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_200Hz', ...  % 4 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_5000Hz', ... % 5 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_500Hz', ...  % 6 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_50Hz', ...   % 7 
        '289x289p_0-350kHz_CHP160_x3_2Vpp_200Hz', ...   % 8
        '289x289p_0-350kHz_CHP20_x3_18Vpp_50Hz', ...    % 9
        '289x289p_0-350kHz_CHP40_x3_18Vpp_50Hz', ...    % 10 
        '289x289p_0-350kHz_CHP80_x3_18Vpp_50Hz', ...    % 11 
        '289x289p_HANN500_x3_18Vpp_50Hz', ...           % 12 
        '289x289p_10HC400kHz_x3_7Vpp_50Hz', ...         % 13
        '289x289p_HANN100_x30_10Vpp_200Hz', ...         % 14 
        '289x289p_HANN100_x3_10Vpp_200Hz', ...          % 15 
        '289x289p_HANN25_x3_10Vpp_200Hz', ...           % 16 
        '289x289p_HANN50_x3_10Vpp_200Hz', ...           % 17 
        '289x289p_HANN50_x3_15Vpp_200Hz', ...           % 18 
        '289x289p_HANN50_x3_20Vpp_200Hz', ...           % 19 
        '289x289p_HANN50_x3_5Vpp_200Hz', ...            % 20 
        '493x493p_0-350kHz_CHP160_x3_18Vpp_50Hz', ...   % 21 
        '493x493p_HANN100_x10_10Vpp_200Hz', ...         % 22 
        '493z493p_0-350kHz_CHP20N_x3', ...              % 23
        '493z493p_0-350kHz_CHP40N_x3', ...              % 24 
        '493z493p_0-500kHz_CHP_40N_x3'};                % 25

disp('Dispersion compensation');
exp_filename = list{exp_test_case};
data_path=fullfile( projectroot, 'data','processed','exp',filesep );
% filename of parameter data
interim_filename_param = ['polar_interim_',exp_filename,'_KXKYF_param']; 
load([data_path,interim_filename_param]); % wavenumber_max fmax beta number_of_wavenumber_points

 % load raw experimental data file
disp('loading data');
load([raw_data_path,exp_filename]); % Data Length Width time
[nY,nX,nft]=size(Data);
nft_padded = 2*nft;
n_padded = nft_padded - nft;
Nexc=n_padded;
Data_padded = zeros(nY,nX,nft_padded);
%Data_padded(1:nY,1:nX,1:nft)=Data; % padding at the end
Data_padded(1:nY,1:nX,n_padded+1:end)=Data;% padding at the begining
dt = time(3)-time(2);
Fs = 1/dt; % Fs - sampling frequency, Units: [Hz]
tt = time(end);  % total calculation time [s]

folder  = model_input_path;
nFile   = length(num_test_case);
success = false(1, nFile);
j=1; % beta=0
for k = num_test_case
    input_name = [model_input_path,num2str(k),'output'];
    load(input_name); % FREQ CG wavenumber
    processed_filename = ['interim_',exp_filename,'_compens']; % filename of processed .mat data
    % check if already exist
    if(overwrite||(~overwrite && ~exist([output_path,processed_filename,'.mat'], 'file')))
        try 
            % load numerical SASE dispersion curves
            input_name = [model_input_path,num2str(k),'output'];
            load(input_name); % FREQ CG wavenumber
            %% PROCESS DATA
            obj_score = 0;
            for j=1:1%length(beta)
                fvec1=squeeze(FREQ(1,:,j)); % mode 1, angle j [Hz]
                fvec2=squeeze(FREQ(2,:,j)); % mode 2, angle j
                fvec3=squeeze(FREQ(3,:,j)); % mode 3, angle j
                fvec4=squeeze(FREQ(4,:,j)); % mode 4, angle j
                kvec=squeeze(wavenumber(:,j)); % angle j [rd/m]
                freq = [Fs/2*linspace(0,1,round(nft_padded)/2)]'; % liearly equally spaced vector of frequencies
                s=squeeze(Data_padded(I,J+K,:)); % signal for dispersion compensation
                figure;
                plot(s);
                y=abs(fft(s));
                y=y(1:nft_padded/2,:);
                figure;
                plot(freq(:,1),y(:,1));
                title('Frequncy components of excitation signal');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Mode 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % interpolate on liearly equally spaced vector of frequencies
                Wavenumber1 = interp1(fvec1,kvec,freq,'linear','extrap'); 
                % post compensation
                sp1(:,1) = designed_waveform(s,-L,freq,Wavenumber1,D0,Nb);
                figure;plot(sp1,'k');
                % draw window
                smax = max(sp1(:,1));
                smin = min(sp1(:,1));
                
                line([(Nexc+1),(Nexc+w),(Nexc+w),(Nexc+1),(Nexc+1)],[smin, smin, smax,smax,smin],'Color','m');
                title('Mode 1');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Mode 2
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % interpolate on liearly equally spaced vector of frequencies
                Wavenumber2 = interp1(fvec2,kvec,freq,'linear','extrap');
                % post compensation
                sp2(:,1) = designed_waveform(s,-L,freq,Wavenumber2,D0,Nb);
                figure;plot(sp2,'k');
                % draw window
                smax = max(sp2(:,1));
                smin = min(sp2(:,1));
         
                line([(Nexc+1),(Nexc+w),(Nexc+w),(Nexc+1),(Nexc+1)],[smin, smin, smax,smax,smin],'Color','m');
                title('Mode 2');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Mode 3
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % interpolate on liearly equally spaced vector of frequencies
                Wavenumber3 = interp1(fvec3,kvec,freq,'linear','extrap');
                % post compensation
                sp3(:,1) = designed_waveform(s,-L,freq,Wavenumber3,D0,Nb);
                figure;plot(sp3,'k');
                % draw window
                smax = max(sp3(:,1));
                smin = min(sp3(:,1));
               
                line([(Nexc+1),(Nexc+w),(Nexc+w),(Nexc+1),(Nexc+1)],[smin, smin, smax,smax,smin],'Color','m');
                title('Mode 3');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Mode 4
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % interpolate on liearly equally spaced vector of frequencies
                Wavenumber4 = interp1(fvec3,kvec,freq,'linear','extrap');
                % post compensation
                sp4(:,1) = designed_waveform(s,-L,freq,Wavenumber3,D0,Nb);
                figure;plot(sp4,'k');
                % draw window
                smax = max(sp4(:,1));
                smin = min(sp4(:,1));
                
                line([(Nexc+1),(Nexc+w),(Nexc+w),(Nexc+1),(Nexc+1)],[smin, smin, smax,smax,smin],'Color','m');
                title('Mode 4');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %sp_sum=sp1+sp2+sp3+sp4;
                sp_sum=abs(sp1)+abs(sp2)+abs(sp3)+abs(sp4);
                figure;plot(sp_sum,'k');
                % draw window
                smax = max(sp_sum(:,1));
                smin = min(sp_sum(:,1));
                
                line([(Nexc+1),(Nexc+w),(Nexc+w),(Nexc+1),(Nexc+1)],[smin, smin, smax,smax,smin],'Color','m');
                title('Sum of modes 1-4');
                obj_score = obj_score + sum(abs(sp_sum(Nexc+1:Nexc+w,1)));
            end
            fprintf('Processing:\n%s\n',filename);
            return;
            %% END OF PROCESSING
            % save data in polar coordinate system
            disp('Saving data...');
            % save processed data to processed data folder
            %save([interim_data_path,processed_filename],'KXKYF_','f_vec','kx_vec', 'ky_vec','-v7.3'); % use compression
            %save([interim_data_path,processed_filename],'KXKYF_','f_vec','kx_vec', 'ky_vec'); % no compression
            [filepath,name,ext] = fileparts(filename);
            fprintf('Successfully processed:\n%s\n', filename);% successfully processed
        catch
            fprintf('Failed: %s\n', exp_filename);
        end
    else
        fprintf('Filename: \n%s \nalready exist\n', processed_filename);
    end
end