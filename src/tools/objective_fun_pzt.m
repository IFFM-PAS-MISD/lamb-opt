function [obj_score] = objective_fun_pzt(signals,L,FREQ,wavenumber,number_of_modes_considered,w,D0,Nb)
% OBJECTIVE_FUN_PZT   score for overlap of numerical and experimental dispersion curves
% 
% Syntax: [obj_score] = objective_fun_pzt(Data_polar,fmax,FREQ,number_of_modes_considered) 
% 
% Inputs: 
%    signals - signals measured by pzt, double, dimensions [number_of_angles,number_of_time_steps]  
%    L - distance from excitation point, double, Units: m 
%    FREQ - Numerical frequency matrix for the same wavenumbers as in experiment, double, 
%           dimensions [number_of_modes,number_of_wavenumber_points,number_of_angles], Units: Hz 
%    wavenumber - vector of wavenumbers, dimensions [1, number_of_wavenumber_points], Units: rad/m
%    number_of_modes_considered - number of modes considered in calculation of the score
% 
% Outputs: 
%    obj_score - Objective function score, double, dimensions [m, n], Units: - 
% 
% Example: 
%    [obj_score] = objective_fun_pzt(Data_polar,fmax,FREQ,number_of_modes_considered)
%    [obj_score] = objective_fun_pzt(Data_polar,fmax,FREQ) 
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

if nargin == 3
    number_of_modes_considered = 4; % default number of modes considered
end
[number_of_angles,nft]=size(signals);
nft_padded = 2*nft;
n_padded = nft_padded - nft;
Nexc=n_padded;
Data_padded = zeros(nX,nft_padded);
Data_padded(:,n_padded+1:end)=Data;% padding at the begining
dt = time(3)-time(2);
Fs = 1/dt; % Fs - sampling frequency, Units: [Hz]  
freq = [Fs/2*linspace(0,1,round(nft_padded)/2)]'; % liearly equally spaced vector of frequencies
           
obj_score = 0;
for j=1:number_of_angles
    kvec=squeeze(wavenumber(:,j)); % angle j [rd/m]
    s=squeeze(Data_padded(j,:)); % signal for dispersion compensation
    for k =1:number_of_modes_considered
        fvec1=squeeze(FREQ(k,:,j)); % mode k, angle j [Hz]

        figure;
        plot(s);
        y=abs(fft(s));
        y=y(1:nft_padded/2,:);
        figure;
        plot(freq(:,1),y(:,1));
        title('Frequncy components of excitation signal');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mode k
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
        title(['Mode ',num2str(k)]);
        obj_score = obj_score + sum(abs(sp1(Nexc+1:Nexc+2,1)));
    end
end

%---------------------- END OF CODE---------------------- 

% ================ [objective_fun_pzt.m] ================  
