function [wavenumber_new] = mode_tracing_new2(FREQ,wavenumber,number_of_modes)
% mode_tracing_new2   Mode tracing using second order approximation
%    sorts roots of Lamb waves dispersion curves (for frequency sweep case) 
% 
% Syntax: [wavenumber_new] = mode_tracing_new2(FREQ,wavenumber,number_of_modes) 
% 
% Inputs: 
%    wavenumber - wavenumbers, double, dimensions
%    [n,number_of_frequency_points], Units: [rad/s] or [Hz]
%    FREQ - frequency vector (equidistant points in frequency
%    domain), double, Units: [Hz] 
%    number_of_modes - number of modes included in sorting, integer
%    (number_of_modes < n)
% 
% Outputs: 
%    wavenumber_new - sorted wavenumbers, double, dimensions
%    [number_of_frequency_points,num_of_modes], Units: [rad/m] or [1/m]
% 
% Example: 
%    [wavenumber_new] = mode_tracing_new2(FREQ,wavenumber,number_of_modes) 
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
FREQ=FREQ/1e3;
[n,number_of_frequency_points]=size(wavenumber);
fmode=zeros(number_of_frequency_points,number_of_modes);
for k=1:number_of_modes
    win=2;
    fmode(1:win,k)=squeeze(wavenumber(k,1:win))';
    for i=1:number_of_frequency_points-win

        finit=zeros(win+1,1);
        kinit=zeros(win+1,1);
        finit(1:win,1)=fmode((i-1)+1:(i-1)+win,k);
        kinit(1:win+1,1)=FREQ((i-1)+1:i+win,1);
        ferror=zeros(number_of_modes,1);
        for j=1:number_of_modes % loop over modes
            finit(win+1,1)=squeeze(wavenumber(j,i+win))';
            
            [p] = polyfit(kinit,finit,2);
            fapprox = polyval(p,kinit);
            ferror(j)=sum((fapprox-finit).^2);
        end
        [A,I]=min(ferror);
        fmode(i+win,k)=squeeze(wavenumber(I,i+win))';
    end
end
wavenumber_new=fmode;
%---------------------- END OF CODE---------------------- 

% ================ [mode_tracing_new2.m] ================  
