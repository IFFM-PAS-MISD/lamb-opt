clear all;close all;
no_of_cycles =3; % best 2 or 2.5 cycles
excit_frequency=250; % [kHz] % 
D0 = 10e3;%  cut off frequency for high-pass Butterworth filter, double , Units: [Hz]
D1 = 300e3;%  cut off frequency for low-pass Butterworth filter, double , Units: [Hz]
Nb=3; % Butterworth filter order, integer

overwrite = false; % allow overwriting existing results if true
number_of_modes_considered = 1; % number of modes considered in calculation of objective function score
number_of_modes=7;
%selected_mode=[2,2,2,3,3,3,3]; % for manual mode tracing
%selected_mode=[3,2,2,3,3,3,3]; % for automatic mode tracing (S0 mode) score 0.9326, 2 cycles
%selected_mode=[1,1,1,1,1,1,1]; % for automatic mode tracing (A0 mode) score 6.5316, 2.5 cycles, 50 kHz
selected_mode=[4,4,4,4,4,4,4]; % for automatic mode tracing (A1 mode)score 0.0086, 3 cycles
%selected_mode=[5,5,5,5,5,5,5]; % for automatic mode tracing (A2 mode) score 0.0282, 3 cycles
L=0.4; % distance between actuator and sensor [m]


%load(['/pkudela_odroid_sensors/lamb_opt/pzt_circ_array_CFRP_uni_Cedrat2/averaged/',num2str(no_of_cycles),'_cycles_',num2str(excit_frequency),'kHz/niscope_avg_waveform.mat']);
%load(['/pkudela_odroid_sensors/lamb_opt/pzt_circ_array_CFRP_uni_Cedrat3/averaged/',num2str(no_of_cycles),'_cycles_',num2str(excit_frequency),'kHz/niscope_avg_waveform.mat']);
load(['/pkudela_odroid_sensors/lamb_opt/pzt_circ_array_CFRP_uni_Cedrat_A0/averaged/',num2str(no_of_cycles),'_cycles_',num2str(excit_frequency),'kHz/niscope_avg_waveform.mat']);
%load(['/pkudela_odroid_sensors/lamb_opt/pzt_circ_array_CFRP_uni_Cedrat_S0/averaged/',num2str(no_of_cycles),'_cycles_',num2str(excit_frequency),'kHz/niscope_avg_waveform.mat']);


load('/home/pkudela/work/projects/opus15/lamb-opt/data/processed/num/genetic_algorithm/ga_unidirectional_C_tensor_known_mass_mut_rnd_offspring_2lay6_out/opt_dispersion_curves.mat');
w=round(1.2*sampleRate/(excit_frequency*1e3/no_of_cycles));% window size in points
dt=1/sampleRate;
%niscope_avg_waveform(1:160,:)=0;
%niscope_avg_waveform(161:end,:)=0;
niscope_avg_waveform(3000:end,1)=0;
niscope_avg_waveform(3000:end,2)=0;
niscope_avg_waveform(3000:end,3)=0;
niscope_avg_waveform(1200:end,4)=0;
niscope_avg_waveform(950:end,5)=0;
%niscope_avg_waveform(3000:end,5)=0;
niscope_avg_waveform(670:end,6)=0;
niscope_avg_waveform(850:end,7)=0;
% for higher frequencies
N=size(niscope_avg_waveform,1);
%niscope_avg_waveform(round(N/2)+1:end,:)=0; % for higher frequencies shorter signal is enough
niscope_avg_waveform(:,7)=niscope_avg_waveform(:,7)/400; % for higher
%frequencies differences in amlitudes at 0 and 90 deg increases
niscope_avg_waveform(:,5)=niscope_avg_waveform(:,5)*2;
n=N;% take whole signal
time=[1:N]*dt-dt;
signals=zeros(7,N);

signals(:,1:N)=niscope_avg_waveform(:,[1,2,3,4,5,6,7])';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% mode tracing automatic
for k=1:7
    [FREQ_new] = mode_tracing_new(squeeze(FREQ(:,:,k)),wavenumber(:,k),number_of_modes);
    FREQ(1:number_of_modes,:,k)=FREQ_new';
end
%figure;plot(squeeze(FREQ(2,:,7)));hold on; plot(squeeze(FREQ(3,:,7)),'r');plot(squeeze(FREQ(4,:,7)),'g');plot(squeeze(FREQ(5,:,7)),'k');
% fmode=zeros(512,number_of_modes);
% for k=2:number_of_modes
%     win=8;
%     fmode(1:win,k)=squeeze(FREQ(k,1:win,1))'/1e3;
%     for i=1:512-win
%         finit=zeros(win+1,1);
%         finit(1:win,1)=fmode((i-1)+1:(i-1)+win,k);
%         kinit(1:win+1,1)=wavenumber((i-1)+1:i+win,1);
%         ferror=zeros(number_of_modes,1);
%         for j=1:number_of_modes % loop over modes
%             finit(win+1,1)=squeeze(FREQ(j,i+win,1))'/1e3;
%             [p] = polyfit(kinit,finit,2);
%             fapprox = polyval(p,kinit);
%             ferror(j)=sum((fapprox-finit).^2);
%         end
%         [A,I]=min(ferror);
%         fmode(i+win,k)=squeeze(FREQ(I,i+win,1))'/1e3;
%     end
% end
% figure;
% plot(fa,'r'); hold on;
% plot(fb,'g');
% plot(fc,'m');
% plot(fd,'c');
% plot(fmode,'k--');

% manual mode tracing
% % angle 1
%     temp1=squeeze(FREQ(2,1:198,1));
%     FREQ(2,1:198,1)=FREQ(3,1:198,1);
%     FREQ(3,1:198,1)=temp1;clear temp1;
%     temp1=FREQ(3,407:end,1);
%     FREQ(3,407:end,1)=FREQ(4,407:end,1);
%     FREQ(4,407:end,1)=temp1;clear temp1;
%     temp1=FREQ(4,1:112,1);
%     FREQ(4,1:112,1)=FREQ(5,1:112,1);
%     FREQ(5,1:112,1)=temp1;clear temp1;
%     temp1=FREQ(5,316:end,1);
%     FREQ(5,316:end,1)=FREQ(4,316:end,1);
%     FREQ(4,316:end,1)=temp1;clear temp1;
% % angle 2
%     temp1=FREQ(3,140:end,2);
%     FREQ(3,140:end,2)=FREQ(4,140:end,2);
%     FREQ(4,140:end,2)=temp1; clear temp1;
% % angle 3
%     temp1=FREQ(4,49:end,3);
%     FREQ(4,49:end,3)=FREQ(3,49:end,3);
%     FREQ(3,49:end,3)=temp1; clear temp1;
% % angle 4
%     temp1=FREQ(4,28:end,4);
%     FREQ(4,28:end,4)=FREQ(3,28:end,4);
%     FREQ(3,28:end,4)=temp1; clear temp1;
% % angle 5
%     temp1=FREQ(4,27:end,5);
%     FREQ(4,27:end,5)=FREQ(3,27:end,5);
%     FREQ(3,27:end,5)=temp1; clear temp1;
% % angle 6
%     temp1=FREQ(4,27:end,6);
%     FREQ(4,27:end,6)=FREQ(3,27:end,6);
%     FREQ(3,27:end,6)=temp1; clear temp1;
%     % angle 7
%     temp=squeeze(FREQ(4,27:end,7));
%     FREQ(4,27:end,7)=FREQ(3,27:end,7);
%     FREQ(3,27:end,7)=temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FREQ must be in angular frequency [rad/s] and wavenumber in [rad/m]
%[score,sig_compens] = objective_fun_pzt2(time(1:n),signals(:,1:n),L,CG,FREQ,wavenumber,number_of_modes_considered,w,D0,Nb);
%[score,sig_compens] = objective_fun_pzt_selected_mode(time(1:n),signals(:,1:n),L,CG,FREQ,wavenumber,selected_mode,w,D0,Nb);
[score,sig_compens] = objective_fun_pzt_selected_mode2(time(1:n),signals(:,1:n),L,CG,FREQ,wavenumber,selected_mode,w,D0,D1,Nb);
score=score*10;
figure;
sp=sum(sig_compens.^2,2);
plot(sp,'k');
smax = max(sp);
smin = min(sp);
line([(N+1),(N+w),(N+w),(N+1),(N+1)],[smin, smin, smax,smax,smin],'Color','m');

%
% frequency plot
fs=sampleRate;
df=fs/N;                            %frequency resolution
sampleIndex = -N/2:N/2-1;   %ordered index for FFT plot
f=sampleIndex*df; 
a=max(abs(fftshift(fft(niscope_avg_waveform(:,5)))));
b=max(abs(fftshift(fft(nifgen_avg_waveform_ch1))));
figure;plot(f,abs(fftshift(fft(niscope_avg_waveform(:,5)))),'r');
xlim([0,1e6]);hold on;
plot(f,a/b*abs(fftshift(fft(nifgen_avg_waveform_ch1))),'b');
