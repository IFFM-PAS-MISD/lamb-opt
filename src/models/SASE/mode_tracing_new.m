function [FREQ_new] = mode_tracing_new(FREQ,wavenumber,number_of_modes)
% mode_tracing_new   Mode tracing using second order approximation
%    sorts roots of Lamb waves dispersion curves  
% 
% Syntax: [FREQ_new] = mode_tracing_new(FREQ,wavenumber,number_of_modes) 
% 
% Inputs: 
%    FREQ - frequencies, double, dimensions
%    [n,number_of_wavenumber_points], Units: [rad/s] or [Hz]
%    wavenumber - wavenumber vector (equidistant points in wavenumber
%    domain), double, Units: [1/m] or [1/rad]
%    number_of_modes - number of modes included in sorting, integer
%    (number_of_modes < n)
% 
% Outputs: 
%    FREQ_new - sorted frequencies, double, dimensions
%    [number_of_wavenumber_points,num_of_modes], Units: [rad/s] or [Hz]
% 
% Example: 
%    [FREQ_new] = mode_tracing_new(FREQ,wavenumber,number_of_modes) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also:  
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 


%% mode tracing automatic
[n,number_of_wavenumber_points]=size(FREQ);
fmode=zeros(number_of_wavenumber_points,number_of_modes);
for k=1:number_of_modes
    win=8;
    fmode(1:win,k)=squeeze(FREQ(k,1:win))'/1e3;
    for i=1:512-win
        finit=zeros(win+1,1);
        finit(1:win,1)=fmode((i-1)+1:(i-1)+win,k);
        kinit(1:win+1,1)=wavenumber((i-1)+1:i+win,1);
        ferror=zeros(number_of_modes,1);
        for j=1:number_of_modes % loop over modes
            finit(win+1,1)=squeeze(FREQ(j,i+win))'/1e3;
            [p] = polyfit(kinit,finit,2);
            fapprox = polyval(p,kinit);
            ferror(j)=sum((fapprox-finit).^2);
        end
        [A,I]=min(ferror);
        fmode(i+win,k)=squeeze(FREQ(I,i+win))'/1e3;
    end
end
FREQ_new=fmode*1e3;
%---------------------- END OF CODE---------------------- 

% ================ [mode_tracing_new.m] ================  
