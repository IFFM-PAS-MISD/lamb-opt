% plot numerical dispersion curves

clear all; close all;

%set(0,'defaulttextinterpreter','none');
%% excitation signals parameters
tt=500/1e6; % total time[s]
nft = 1024*4; % number of sampling points
dt = tt/nft;
% narrow_band
nc_narrow = 8; % number of cycles
f_1_narrow = 200e3/nc_narrow; % modulation frequency [Hz]
f_2_narrow = nc_narrow*f_1_narrow; % carrier frequency [Hz]
t_1_narrow = 10/1e6; % initiation time [s]
w_narrow = round(1/dt/f_1_narrow);% window size
% wide_band
nc_wide = 1.5; % number of cycles
f_1_wide = 20e3/nc_wide; % modulation frequency [Hz]
f_2_wide = nc_wide*f_1_wide; % carrier frequency [Hz]
t_1_wide = 10/1e6; % initiation time [s]
w_wide = round(1/dt/f_1_wide);% window size
%%
% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true

% figure parameters
% size 12cm by 8cm (1-column text)
%fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
fig_width = 7; fig_height = 5; 
% create path to the numerical model data folder
modelfolder = 'SASE';
modelname = 'SASE2';
modelname2 = 'dispersion_effect';
radians = false;
% create output path
output_path = prepare_figure_paths(modelfolder,modelname2);

% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname,'_out'], filesep );
% create path to the experimental processed data folder
exp_input_data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
% and conversion to polar coordinate system
filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_param'; 

% load experimental parameter data file
load([exp_input_data_path,filename]); % wavenumber_max fmax beta number_of_wavenumber_points
number_of_frequency_points = number_of_wavenumber_points;
number_of_angles = length(beta);
fvec = linspace(0,fmax,number_of_frequency_points);
if(~radians)
   wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
end
fprintf('Making figures: %s\n', modelname2);
%% plot excitation signals
fr=zeros(nft,1);
for n=1:nft
  fr(n)=(n-1)/tt; % frequency [Hz]
end
% narrow band
figfilename='excitation_narrow_time';
[t_narrow,st_narrow]=Hanning_signal(dt,nft,f_1_narrow,f_2_narrow,t_1_narrow);
    figure;
    plot(t_narrow*1e6,st_narrow,'linewidth',0.5);
    box on;
    set(gcf,'Color','white');
    %axis([0 60 -1.05 1.05]);
    axis([0 (2*t_1_narrow+t_narrow(w_narrow))*1e6 -1.05 1.05]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$t$ [$\mu$s]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$A$'},'Fontsize',12,'interpreter','latex'); 
    set(gca,'FontName','Times');
    fig = gcf; 
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 
Pyy_narrow = fft(st_narrow,nft)*2/nft;
Pyy_narrow_norm=abs(Pyy_narrow)/max(abs(Pyy_narrow));
figfilename='excitation_narrow_frequency';
    figure;
    plot(fr/1e3,abs(Pyy_narrow_norm),'linewidth',0.5,'color',[0.850,0.325,0.098]);
    box on;
    set(gcf,'Color','white');
    %axis([0 400 0 1.05]);
    axis([0 300 0 1.05]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$A$'},'Fontsize',12,'interpreter','latex'); 
    set(gca,'FontName','Times');
    fig = gcf; 
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 
% wide band
figfilename='excitation_wide_time';
[t_wide,st_wide]=Hanning_signal(dt,nft,f_1_wide,f_2_wide,t_1_wide);
    figure;
    plot(t_wide*1e6,st_wide,'linewidth',0.5);
    box on;
    set(gcf,'Color','white');
    %axis([0 40 -1.05 1.05]);
    axis([0 (2*t_1_wide+t_wide(w_wide))*1e6 -1.05 1.05]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$t$ [$\mu$s]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$A$'},'Fontsize',12,'interpreter','latex'); 
    set(gca,'FontName','Times');
    fig = gcf; 
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 
Pyy_wide = fft(st_wide,nft)*2/nft;
Pyy_wide_norm=abs(Pyy_wide)/max(abs(Pyy_wide));
figfilename='excitation_wide_frequency';
    figure;
    plot(fr/1e3,abs(Pyy_wide_norm),'linewidth',0.5,'color',[0.850,0.325,0.098]);
    box on;
    set(gcf,'Color','white');
    axis([0 300 0 1.05]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$A$'},'Fontsize',12,'interpreter','latex'); 
    set(gca,'FontName','Times');
    fig = gcf; 
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 
    
%% plot dispersion curves with overlapped frequency range
input_name = [model_input_path,num2str(7),'output'];
load(input_name); % FREQ CG wavenumber
FREQ_A0 = FREQ(1,:,3);
CG_A0 = CG(1,:,3);
if(~radians)
    wavenumber = wavenumber/(2*pi); % linear scale [1/m]
end
kvec=squeeze(wavenumber(:,3)); % angle j
% size 7cm by 5cm (2-column slides)
fig_width = 7; fig_height = 5;
X=[160,240,240,160,160];
Y=[0,0,1600,1600,0];
% narrow band
figure;
figfilename = 'A0_dispersion_less_dispersive';          
    fill(X,Y,[0.8,0.8,0.8]);
    hold on;
    yyaxis left;
    plot(FREQ_A0(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[0,0,0]);
    hold on;
    box on;
    set(gcf,'Color','white');
    axis([0 350 0 min(wavenumber_max)]);
    %axis([0 500 0 min(wavenumber_max)]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
    if(radians)
        ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
    else
        ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
    end
    yyaxis right;
    %plot(fr/1e3,abs(Pyy_narrow_norm)*250,'linewidth',0.5,'color',[0.635,0.078,0.184]);
    plot(fr/1e3,abs(Pyy_narrow_norm),'linewidth',1,'color',[0.850,0.325,0.098]);
    ylim([0 1]);
    set(gca,'FontName','Times');
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 
figure;
figfilename = 'A0_phase_velocity_less_dispersive';          
    fill(X,Y,[0.8,0.8,0.8]);
    hold on;
    yyaxis left;
    plot(FREQ_A0(2:end)/1e3,FREQ_A0(2:end)./kvec(2:end)','linewidth',1,'color',[0,0,0]);
    hold on;
    box on;
    set(gcf,'Color','white');
    axis([0 350 0 1500]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$c_p$ [m/s]'},'Fontsize',12,'interpreter','latex');
    yyaxis right;
    %plot(fr/1e3,abs(Pyy_narrow_norm)*250,'linewidth',0.5,'color',[0.635,0.078,0.184]);
    plot(fr/1e3,abs(Pyy_narrow_norm),'linewidth',1,'color',[0.850,0.325,0.098]);
    ylim([0 1]);
    set(gca,'FontName','Times');
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 

figure;
figfilename = 'A0_group_velocity_less_dispersive';          
    fill(X,Y,[0.8,0.8,0.8]);
    hold on;
    yyaxis left;
    plot(FREQ_A0(2:end)/1e3,CG_A0(2:end),'linewidth',1,'color',[0,0,0]);
    hold on;
    box on;
    set(gcf,'Color','white');
    axis([0 350 0 1600]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$c_g$ [m/s]'},'Fontsize',12,'interpreter','latex');
    yyaxis right;
    %plot(fr/1e3,abs(Pyy_narrow_norm)*250,'linewidth',0.5,'color',[0.635,0.078,0.184]);
    plot(fr/1e3,abs(Pyy_narrow_norm),'linewidth',1,'color',[0.850,0.325,0.098]);
    ylim([0 1]);
    set(gca,'FontName','Times');
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 
% wide band
X=[0,42,42,0,0];
figure;
figfilename = 'A0_dispersion_dispersive';          
    fill(X,Y,[0.8,0.8,0.8]);
    hold on;
    yyaxis left;
    plot(FREQ_A0(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[0,0,0]);
    hold on;
    box on;
    set(gcf,'Color','white');
    axis([0 350 0 min(wavenumber_max)]);
    %axis([0 500 0 min(wavenumber_max)]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
    if(radians)
        ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
    else
        ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
    end
    yyaxis right;
    %plot(fr/1e3,abs(Pyy_wide_norm)*250,'linewidth',0.5,'color',[0.635,0.078,0.184]);
    plot(fr/1e3,abs(Pyy_wide_norm),'linewidth',1,'color',[0.850,0.325,0.098]);
    set(gca,'FontName','Times');
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600');   
figure;
figfilename = 'A0_phase_velocity_dispersive';          
    fill(X,Y,[0.8,0.8,0.8]);
    hold on;
    yyaxis left;
    plot(FREQ_A0(2:end)/1e3,FREQ_A0(2:end)./kvec(2:end)','linewidth',1,'color',[0,0,0]);
    hold on;
    box on;
    set(gcf,'Color','white');
    axis([0 350 0 1500]);
    %axis([0 500 0 min(wavenumber_max)]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$c_p$ [m/s]'},'Fontsize',12,'interpreter','latex');
    yyaxis right;
    %plot(fr/1e3,abs(Pyy_wide_norm)*250,'linewidth',0.5,'color',[0.635,0.078,0.184]);
    plot(fr/1e3,abs(Pyy_wide_norm),'linewidth',1,'color',[0.850,0.325,0.098]);
    set(gca,'FontName','Times');
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600');   
figure;
figfilename = 'A0_group_velocity_dispersive';
    fill(X,Y,[0.8,0.8,0.8]);
    hold on;
    yyaxis left;
    plot(FREQ_A0(2:end)/1e3,CG_A0(2:end),'linewidth',1,'color',[0,0,0]);
    hold on;
    box on;
    set(gcf,'Color','white');
    axis([0 350 0 1600]);
    %axis([0 500 0 min(wavenumber_max)]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$c_g$ [m/s]'},'Fontsize',12,'interpreter','latex');
    yyaxis right;
    %plot(fr/1e3,abs(Pyy_wide_norm)*250,'linewidth',0.5,'color',[0.635,0.078,0.184]);
    plot(fr/1e3,abs(Pyy_wide_norm),'linewidth',1,'color',[0.850,0.325,0.098]);
    set(gca,'FontName','Times');
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600');   
 %% dispersion effect - narrow band
L=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]*4; % distance of compensation
nft=50000;   % number of points in the signal
%nft=20000;   % number of points in the signal
fm=f_1_narrow;    % fm - modulation frequency, double, Units: [Hz] 
fc=f_2_narrow;      % frequency of the carrier signal [Hz]
Nexc = nft/2;%10000; % 
Ap=1; % peak amplitude
D0 = 10e3;%  cut off frequency for high-pass Butterworth filter, double , Units: [Hz]
Nb=1; % Butterworth filter order, integer
%%
%
% min sampling 6 MHz
% f_2 = 70e3; sampling =  f_2*20; >15
% f_2 = 100e3; sampling =  f_2*15; >15
% f_2 = 500e3; sampling =  f_2*15; >12

Fs = fc*15; % Fs - sampling frequency, Units: [Hz]
dt = 1/Fs;
tt = (nft-1)*dt;  % total calculation time [s]
t_1= nft/2*dt;       % excitation initiation time (delay) [s]

w = round(Fs/fm);% window size
% hanning windowed signal
tic;
[t,st]=Hanning_sig(dt,nft,fm,fc,t_1,w);
st=-Ap*st;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode 1
Frequency=FREQ_A0'; % Frequency components of dispersion curve [Hz]
f = Fs/2*linspace(0,1,round(nft/2)); % liearly equally spaced vector of frequencies
Wavenumber = interp1(Frequency,kvec,f','linear','extrap');
[I]=find(f<20);
f(I)=0;
Wavenumber(I)=0;
[I2]=find(Wavenumber<0);
f(I2)=0;
Wavenumber(I2)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L1=L;
M = length(L); % number of signals
st = repmat(st,[1,M]);
f = repmat(f',[1,M]);
L= repmat(L,[nft/2,1]);
Wavenumber = repmat(Wavenumber,[1,M]);
% forward propagating of excitation signal d=L
d=L; % resembles acquired signal at pzt (reflected from damage)

mode1_L = designed_waveform(st,d,f,Wavenumber,D0,Nb);
fig_width = 12; fig_height = 4; 
for k=1:M
    figfilename=['dispersion_effect_less_dispersive_L_',num2str((k))];
    figure;
    plot((t-t_1)*1e6,mode1_L(:,k),'linewidth',0.5);
    box on;
    set(gcf,'Color','white');
    %axis([0 200 -1.05 1.05]);
    axis([0 500 -1.05 1.05]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$t$ [$\mu$s]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$A$'},'Fontsize',12,'interpreter','latex'); 
    title({['Propagation distance: ',num2str(L1(k)),' m']},'Fontsize',12,'interpreter','latex');
    set(gca,'FontName','Times');
    fig = gcf; 
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 
end
%% dispersion effect - wide band
L=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]*4; % distance of compensation
nft=50000;   % number of points in the signal
%nft=20000;   % number of points in the signal
fm=f_1_wide;    % fm - modulation frequency, double, Units: [Hz] 
fc=f_2_wide;      % frequency of the carrier signal [Hz]
Nexc = nft/2;%10000; % 
Ap=1; % peak amplitude
D0 = 1e3;%  cut off frequency for high-pass Butterworth filter, double , Units: [Hz]
Nb=1; % Butterworth filter order, integer
%%
%
% min sampling 6 MHz
% f_2 = 70e3; sampling =  f_2*20; >15
% f_2 = 100e3; sampling =  f_2*15; >15
% f_2 = 500e3; sampling =  f_2*15; >12

Fs = fc*15; % Fs - sampling frequency, Units: [Hz]
dt = 1/Fs;
tt = (nft-1)*dt;  % total calculation time [s]
t_1= nft/2*dt;       % excitation initiation time (delay) [s]

w = round(Fs/fm);% window size
% hanning windowed signal
tic;
[t,st]=Hanning_sig(dt,nft,fm,fc,t_1,w);
st=-Ap*st;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode 1
Frequency=FREQ_A0'; % Frequency components of dispersion curve [Hz]
f = Fs/2*linspace(0,1,round(nft/2)); % liearly equally spaced vector of frequencies
Wavenumber = interp1(Frequency,kvec,f','linear','extrap');
[I]=find(f<20);
f(I)=0;
Wavenumber(I)=0;
[I2]=find(Wavenumber<0);
f(I2)=0;
Wavenumber(I2)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L1=L;
M = length(L); % number of signals
st = repmat(st,[1,M]);
f = repmat(f',[1,M]);
L= repmat(L,[nft/2,1]);
Wavenumber = repmat(Wavenumber,[1,M]);
% forward propagating of excitation signal d=L
d=L; % resembles acquired signal at pzt (reflected from damage)

mode1_L = designed_waveform(st,d,f,Wavenumber,D0,Nb);
fig_width = 12; fig_height = 4; 
for k=1:M
    figfilename=['dispersion_effect_dispersive_L_',num2str((k))];
    figure;
    plot((t-t_1)*1e6,mode1_L(:,k),'linewidth',0.5);
    box on;
    set(gcf,'Color','white');
    %axis([0 200 -1.05 1.05]);
    axis([0 800 -1.05 1.05]);
    set(gca, 'Layer', 'Top');
    set(gca,'Fontsize',10,'linewidth',1);
    xlabel({'$t$ [$\mu$s]'},'Fontsize',12,'interpreter','latex');
    ylabel({'$A$'},'Fontsize',12,'interpreter','latex'); 
     title({['Propagation distance: ',num2str(L1(k)),' m']},'Fontsize',12,'interpreter','latex');
    set(gca,'FontName','Times');
    fig = gcf; 
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    print([output_path,figfilename],'-dpng', '-r600'); 
end
%% END PLOTTING

