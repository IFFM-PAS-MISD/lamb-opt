function [stt] = designed_waveform(st,L,f,Wavenumber,D0,Nb)
% DESIGNED_WAVEFORM   Apply predispersion to original signal 
% 
% Syntax: [stt] = designed_waveform(st,L,om,wavenumber,D0,Nb) 
% 
% Inputs: 
%    st - time domain signal, matrix of doubles, dimensions [NFFT,M]
%         NFFT=length(st),M - number of singnals 
%    L - compensation distance, vector of doubles, dimensions [1,M], Units: [m]
%    f - frequency components, repeated M colums, dimensions [NFFT,M], Units: [Hz]
%    Wavenumber - wavenumber compponents of dispersion curve wavenumber(f),
%                 repeated M columns, double, dimensions [NFFT,M], Units: rad/s
%    D0 - cut off frequency for high-pass Butterworth filter, double , Units: [Hz]
%    Nb - Butterworth filter order, integer
% 
% Outputs: 
%    stt - precompensated signal (i.e. designed excitation waveform)
%          matrix of doubles, dimensions [NFFT,M]
% 
% Example: 
%    [stt] = designed_waveform(st,L,f,Wavenumber,D0,Nb)
%    [stt] = designed_waveform(st,0.5,f,Wavenumber,8e3,1)
%    [stt] = designed_waveform(st,1.0,f,Wavenumber,10000,2) 
% 
% Other m-files required:  butterworth_high_pass1D.m 
% MAT-files required: none 
% See also:  BUTTERWORTH_HIGH_PASS1D 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

NFFT=length(st);
% fft of the signal
Y = fft(st); %

Y = Y(1:round(NFFT/2),:); % select part up to Nyquist frequency

% apply high pass filter >10kHz 
[H] = butterworth_high_pass1D(f,D0,Nb);

Y=H.*Y;

% apply precompensation
G = Y.*exp(-1i.*Wavenumber.*L);

% Method 1 (Zhibo) 
% Gr = [real(G(1:end));0;flipud(real(G(2:end)))];
% Gi = [-imag(G(1:end));0;flipud(imag(G(2:end)))];
% Q = Gr-1i*Gi; stt=ifft(Q); % should be Q = Gr-1i*Gi; not Q = Gr+1i*Gi;

% Method 2
%Q=[G;0;conj(flipud(G(2:end)))]; stt=ifft(Q,NFFT);

% Method 3
stt=ifft(G,NFFT,'symmetric');


%---------------------- END OF CODE---------------------- 

% ================ [designed_waveform.m] ================  
