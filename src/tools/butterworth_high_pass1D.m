function [H] = butterworth_high_pass1D(f,D0,Nb)
% BUTTERWORTH_HIGH_PASS1D   Apply Butterworth high pass filter in frequency domain 
%     Frequency componnents corresponds to S - signal in frequency domain S=fft(s) of dimensions [m, 1] 
%
% Syntax: [H] = butterworth_high_pass1D(s,f,D0,n)
% 
% Inputs: 
%    
%    f - frequency components, vector of doubles, dimensions [m, 1],  Units: [Hz]
%    D0 - cut off frequency, double , Units: [Hz]
%    Nb - filter order, integer
% 
% Outputs: 
%    H - Filter mask in frequency domain, vector of doubles, dimensions [m, 1], m = length(f)
% 
% Example: 
%    [H] = butterworth_high_pass1D(S,f,D0,Nb)
%    [H] = butterworth_high_pass1D(S,f,10000,2) 
% 
% Other m-files required: none 
% MAT-files required: none 
% See also: APPLY_REJECT_HIGH_PASS1D 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

H = 1./(1 + (D0./f).^(2*Nb) );

%---------------------- END OF CODE---------------------- 

% ================ [butterworth_high_pass1D.m] ================  
