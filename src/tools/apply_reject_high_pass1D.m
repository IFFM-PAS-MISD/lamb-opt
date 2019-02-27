function [Smod] = apply_reject_high_pass1D(S,f,D0)
% APPLY_REJECT_HIGH_PASS1D   Apply rejecting high pass filter in frequency domain 
% 
% Syntax: [S] = apply_reject_high_pass1D(S,f,D0)
% 
% Inputs: 
%    S - signal in frequency domain S=fft(s) of dimensions [m, 1], m = length(S)
%    f - frequency components, vector of doubles, dimensions [m, 1],  Units: [Hz]
%    D0 - cut off frequency, double , Units: [Hz]
% 
% Outputs: 
%    
%     Smod - modified signal in frequency domain with applied filter, dimensions [m, 1]
% 
% Example: 
%    [S] = apply_reject_high_pass1D(S,f,D0) 
%    [S] = apply_reject_high_pass1D(S,f,10000) 
% 
% Other m-files required: none 
% MAT-files required: none 
% See also: BUTTERWORTH_HIGH_PASS1D 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

I = find(f<D0);
S(I)=0;
Smod=S;
% figure;
% plot(f', abs(y))

%---------------------- END OF CODE---------------------- 

% ================ [apply_reject_high_pass1D.m] ================  
