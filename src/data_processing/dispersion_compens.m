clear all;close all;
no_of_cycles =3; % best 2 or 2.5 cycles
excit_frequency=190; % [kHz] % best 40-60 kHz
D0 = 10e3;%  cut off frequency for high-pass Butterworth filter, double , Units: [Hz]
Nb=1; % Butterworth filter order, integer

overwrite = false; % allow overwriting existing results if true
number_of_modes_considered = 1; % number of modes considered in calculation of objective function score
L=0.4; % distance between actuator and sensor [m]


%load(['/pkudela_odroid_sensors/lamb_opt/pzt_circ_array_CFRP_uni_Cedrat2/averaged/',num2str(no_of_cycles),'_cycles_',num2str(excit_frequency),'kHz/niscope_avg_waveform.mat']);
%load(['/pkudela_odroid_sensors/lamb_opt/pzt_circ_array_CFRP_uni_Cedrat3/averaged/',num2str(no_of_cycles),'_cycles_',num2str(excit_frequency),'kHz/niscope_avg_waveform.mat']);
load(['/pkudela_odroid_sensors/lamb_opt/pzt_circ_array_CFRP_uni_Cedrat_A0/averaged/',num2str(no_of_cycles),'_cycles_',num2str(excit_frequency),'kHz/niscope_avg_waveform.mat']);

load('/home/pkudela/work/projects/opus15/lamb-opt/data/processed/num/genetic_algorithm/ga_unidirectional_C_tensor_known_mass_mut_rnd_offspring_2lay6_out/opt_dispersion_curves.mat');
w=round(1.2*sampleRate/(excit_frequency*1e3/no_of_cycles));% window size in points
dt=1/sampleRate;

% for higher frequencies
N=size(niscope_avg_waveform,1);
%niscope_avg_waveform(round(N/2)+1:end,:)=0; % for higher frequencies shorter signal is enough
%niscope_avg_waveform(:,7)=niscope_avg_waveform(:,7)/10; % for higher
%frequencies differences in amlitudes at 0 and 90 deg increases
n=N;% take whole signal

time=[1:N]*dt-dt;
signals=zeros(7,N);

signals(:,1:N)=niscope_avg_waveform(:,[1,2,3,4,5,6,7])';
% signals(6,:)=signals(6,:)/2;
% signals(7,:)=signals(7,:)/2;
% FREQ must be in angular frequency [rad/s] and wavenumber in [rad/m]
[score,sig_compens] = objective_fun_pzt2(time(1:n),signals(:,1:n),L,CG,FREQ,wavenumber,number_of_modes_considered,w,D0,Nb);
figure;
sp=sum(sig_compens.^2,2);
plot(sp,'k');
smax = max(sp);
smin = min(sp);
line([(N+1),(N+w),(N+w),(N+1),(N+1)],[smin, smin, smax,smax,smin],'Color','m');

fs=sampleRate;
df=fs/N;                            %frequency resolution
sampleIndex = -N/2:N/2-1;   %ordered index for FFT plot
f=sampleIndex*df; 
figure;plot(f,abs(fftshift(fft(niscope_avg_waveform(:,3)))),'r');
xlim([0,1e6]);