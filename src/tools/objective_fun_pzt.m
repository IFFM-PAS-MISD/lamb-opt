function [obj_score] = objective_fun_pzt(time,signals,L,CG,FREQ,wavenumber,number_of_modes_considered,w,D0,Nb)
% OBJECTIVE_FUN_PZT   score for dispersion compensated signals
%    dispersion curve is calculated numerically whearas signals are experimental
% 
% Syntax: [obj_score] = objective_fun_pzt(Data_polar,fmax,FREQ,number_of_modes_considered) 
% 
% Inputs: 
%    time - time vector, double, dimension [number_of_time_steps,1]
%    signals - signals measured by pzt, double, dimensions [number_of_angles,number_of_time_steps]  
%    L - distance from excitation point, double, Units: m 
%    FREQ - Numerical frequency matrix for the same wavenumbers as in experiment, double, 
%           dimensions [number_of_modes,number_of_wavenumber_points,number_of_angles], Units: Hz 
%    wavenumber - vector of wavenumbers, dimensions [1, number_of_wavenumber_points], Units: rad/m
%    number_of_modes_considered - number of modes considered in calculation of the score
%    w - window size in points, integer
%    D0 - cut off frequency for high-pass Butterworth filter, double , Units: [Hz]
%    Nb - Butterworth filter order, integer
% 
% Outputs: 
%    obj_score - Objective function score, double, dimensions [m, n], Units: - 
% 
% Example: 
%    [obj_score] = objective_fun_pzt(time,signals,L,FREQ,wavenumber,number_of_modes_considered,w,D0,Nb)
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: OBJECTIVE_FUN 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 


[number_of_angles,nft]=size(signals);
nft_padded = 2*nft;
%nft_padded = nft;
n_padded = nft_padded - nft;
Nexc=n_padded;
Data_padded = zeros(number_of_angles,nft_padded);
Data_padded(:,n_padded+1:end)=signals;% padding at the begining
dt = time(3)-time(2);
Fs = 1/dt; % Fs - sampling frequency, Units: [Hz]  
freq = [Fs/2*linspace(0,1,round(nft_padded)/2)]'; % liearly equally spaced vector of frequencies
obj_score = 0;
for j=1:number_of_angles
    kvec=squeeze(wavenumber(:,j)); % angle j [rd/m]
    s=squeeze(Data_padded(j,:))'; % signal for dispersion compensation
    sp1 = zeros(length(s),1);
    sp = zeros(length(s),1);
    %wavenumber_step = kvec(2) - kvec(1);
    % mode tracing is unnecessary 
    %[~,FREQ_new] = mode_tracing_pade(CG(1:number_of_modes_considered,:,j),squeeze(FREQ(1:number_of_modes_considered,:,j)),wavenumber_step);
    for k =1:number_of_modes_considered
        fvec1=squeeze(FREQ(k,:,j))'; % mode k, angle j [Hz]
        %fvec1=FREQ_new(k,:)'; % mode k, angle j [Hz]
%         figure;
%         plot(s);
%         y=abs(fft(s));
%         y=y(1:nft_padded/2,:);
%         figure;
%         plot(freq(:,1),y(:,1));
%         title('Frequncy components of excitation signal');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mode k
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % interpolate on liearly equally spaced vector of frequencies
        Wavenumber1 = interp1(fvec1,kvec,freq,'linear','extrap'); 
        [I]=find(freq<20);
        freq(I)=0;
        Wavenumber1(I)=0;
        [I2]=find(Wavenumber1<0);
        freq(I2)=0;
        Wavenumber1(I2)=0;
        % post compensation
        sp(:,1) = designed_waveform(s,-L,freq,Wavenumber1,D0,Nb);
        sp1 = sp1 + sp;
        %figure;plot(sp1,'k');
        % draw window
%         smax = max(sp1(:,1));
%         smin = min(sp1(:,1));
%         line([(Nexc+1),(Nexc+w),(Nexc+w),(Nexc+1),(Nexc+1)],[smin, smin, smax,smax,smin],'Color','m');
%         title(['Mode ',num2str(k)]);      
    end
    obj_score = obj_score + sum(abs(sp1(Nexc+1:Nexc+w,1)));
end

%---------------------- END OF CODE---------------------- 

% ================ [objective_fun_pzt.m] ================  
