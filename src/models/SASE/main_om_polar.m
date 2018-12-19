%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                          main                       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calls SASE function for a given frequency (in kHz) and gets 
% the real and imaginary wave numbers.
% And plots attenuation at different angles (beta) for given eigen mode 
% (eig_no)
% freq in kHz                       beta in degrees
% np = order of elements (3<=np<=5)
% nele_layer = no. of elements per ply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;
% clc;

%% Input for SAFE

% Mode, frequency and angles and order of sprectal elements
beta = 0:5:90; % angles for dispersion curves in polar plot for fixed frequency
np = 3;
nele_layer = 1;
freq=100; % frequency for polar plots [kHz]
theta=0; % angle for dispersion curves
%fmin=1; %minimal frequency [kHz]
wavenumber_min = 1; % minimal wavenumber [1/m]
frmax=500; % maximal frequency for dispersion curves [kHz]
wavenumber_max = 2500; % maximal wavenumber for dispersion curves [1/m]
%number_of_freq_points=500;
number_of_wavenumber_points=1000;
%fr_step=(frmax-fmin)/(number_of_freq_points-1); % frequency step [kHz]
k_step=(wavenumber_max-wavenumber_min)/(number_of_wavenumber_points-1); % wavenumber step [1/m]
nbeta = length(beta);
num_of_modes=9; % number of modes to display on dispersion curves
type='composite'; % type='composite' or type='isotropic'

fname='test';
sorting = 'yes';
%% SAFE
disp('SASE...');

%
%% Input
SASEinput

%% Transform material properties

C = cell(nlayers,1);

for ii=1:nlayers
    theta = rot_angles(ii);
    %C_tmp = transprop(C0r+1i*C0i,stack_dir,1*theta);
    C_tmp = transprop(C0r,stack_dir,1*theta); % no damping, real elastic constants
    C{ii} = C_tmp;
end


%% Get K and M
% Get global stiffness and mass matrices (independent of beta and freq.)
[K, M] = get_KM(nele_layer,np,h,C,rho);
%% dispersion curves
FREQ = zeros(length(M),number_of_wavenumber_points,length(beta));
CG = zeros(length(M),number_of_wavenumber_points,length(beta));
for j=1:length(beta)
    
num_of_modes=length(M);
om_real = zeros(length(M),number_of_wavenumber_points);
om_imag = zeros(length(M),number_of_wavenumber_points);
cg = zeros(length(M),number_of_wavenumber_points);
om_real2 = zeros(length(M),number_of_wavenumber_points);
om_imag2 = zeros(length(M),number_of_wavenumber_points);
cg2 = zeros(length(M),number_of_wavenumber_points);
wavenumber = zeros(1,number_of_wavenumber_points);


% in Np/m -- convert to dB/m -- 1 Np/m = 20/log(10) dB/m
disp('SASE...');
for k=1:number_of_wavenumber_points
    [k number_of_wavenumber_points]
    wavenumber(k) = wavenumber_min+(k-1)*k_step;
    [cg(:,k),mode_shapes, om_real(:,k), om_imag(:,k)] = SASE_om(K,M,beta(j),wavenumber(k));
end


%%
c='bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';
% figure(12); hold on;
% for k=1:num_of_modes
%     s=[c(k),'.'];
%     plot(om_real(k,:)/2/pi/1e3,wavenumber,s,'LineWidth',2);
% end
% Hxl=xlabel('Frequency [kHz]');
% %Hxl=xlabel('Frequency [Hz]');
% Hyl=ylabel('Wavenumber [1/m]');
% 
% set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
% set(gcf,'Color','white');set(gca,'FontSize',11);
% title([num2str(beta(1)),' [deg]']);
% axis([0 500000/1e3 0 1800]);
% 
% 
% figure(13); hold on;
% for k=1:num_of_modes
%     s=[c(k),'.'];
%     plot(om_real(k,:)/2/pi/1e3,om_real(k,:)./wavenumber,s,'LineWidth',2);
% end
% Hxl=xlabel('Frequency [kHz]');
% Hyl=ylabel('Phase velocity');
% set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
% set(gcf,'Color','white');set(gca,'FontSize',11);
% title([num2str(beta(1)),' [deg]']);
% axis([0 500000/1e3 0 10000]);
% 
% figure(14); hold on;
% for k=1:num_of_modes
%     s=[c(k),'.'];
%     plot(om_real(k,:)/2/pi/1e3,(real(cg(k,:))),s,'LineWidth',2);
% end
% Hxl=xlabel('Frequency [kHz]');
% Hyl=ylabel('Group velocity');
% set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
% set(gcf,'Color','white');set(gca,'FontSize',11);
% title([num2str(beta(1)),' [deg]']);
% axis([0 500000/1e3 0 10000]);


%% mode-tracing
% Pade expansion method

om_new = om_real;
cg_new=cg;

for i=1:num_of_modes
    I=i;cc=0;
for k=number_of_wavenumber_points:-1:2
   
    if (cc==0) % first order Taylor approximation
        cc=cc+1;
        cg0=cg(I,k);
        cg1=cg(I,k-1);
        gamma0=om_real(I,k); 
        gamma1=cg0;
        om_t=gamma0-gamma1*k_step; % first order Taylor approximation (backward)
        %cg_t=cg1; 
        cg_t=cg0; 
        delta_om = abs((om_t-om_real(:,k-1))/om_t); % deviation of angular frequency in Taylor method
        delta_cg = abs((cg_t-cg(:,k-1))/cg_t);  % deviation of group velocity in Taylor method
        delta = delta_om.^2+delta_cg.^2;
        [a,I] = min(delta_om);
        %ind(i,k-1) = I;
        om_new(i,k-1) = om_real(I,k-1);
        cg_new(i,k-1) = cg(I,k-1);
    elseif (cc==1) % second order Taylor approximation
        cc=cc+1; 
        cg0=cg(I,k+1);
        cg1=cg(I,k);
        gamma0=om_real(I,k); 
        gamma1=cg0;
        gamma2=1/2*(cg1-cg0)/k_step;
        om_t=gamma0-gamma1*k_step-gamma2*k_step^2; % second order Taylor approximation (backward)
        cg_t=2*cg1-cg0; % linear Taylor approximation
        delta_om = abs((om_t-om_real(:,k-1))/om_t); % deviation of angular frequency in Taylor method
        delta_cg = abs((cg_t-cg(:,k-1))/cg_t);  % deviation of group velocity in Taylor method
        delta = delta_om.^2+delta_cg.^2;
        [a,I] = min(delta_om);
        %ind(i,k-1) = I;
        om_new(i,k-1) = om_real(I,k-1);
        cg_new(i,k-1) = cg(I,k-1);
    else
        cc=cc+1; 
        cg0=cg(I,k+2);
        cg1=cg(I,k+1);
        cg2=cg(I,k);
    
        gamma0=om_real(I,k);
        gamma1=cg2;
        gamma2=1/2*(cg2-cg1)/k_step;
        gamma3=1/6*(cg2-2*cg1+cg0)/(k_step^2);
    
        beta2=(gamma2^2-gamma3*gamma1)/(gamma1^2-gamma2*gamma0);
        beta1=-(gamma2+gamma0*beta2)/gamma1;
        alpha0=gamma0;
        alpha1=gamma1+gamma0*beta1;
    
        %om_p=(alpha0 + alpha1 * k_step)/(1 + beta1*k_step + beta2*k_step^2); %Pade expansion of order [1/2](forward)
        %om_p=(alpha0 - alpha1 * k_step)/(1 -( beta1*k_step + beta2*k_step^2)); % Pade expansion of order [1/2](backward)
        om_t=gamma0-gamma1*k_step-gamma2*k_step^2-gamma3*k_step^3; % Third order Taylor approximation (backward)
    
        cg_t=2*cg2-cg1; % linear Taylor approximation
        
        delta_om = abs((om_t-om_real(:,k-1))/om_t); % deviation of angular frequency in Taylor method
        delta_cg = abs((cg_t-cg(:,k-1))/cg_t);  % deviation of group velocity in Taylor method
        delta = delta_om.^2+delta_cg.^2;
        %[a,I] = min(delta);
        [a,I] = min(delta_om);
        
        %ind(i,k-1) = I;
        om_new(i,k-1) = om_real(I,k-1);
        cg_new(i,k-1) = cg(I,k-1);
        if(I~=i)
            cc=0;
        end
    end

 
end
FREQ(:,:,j) = om_new/2/pi;
CG(:,:,j) = cg_new;

end
%%
end % end of loop over angle beta

%% figures
figure(16); hold on;
for k=1:num_of_modes
    s=[c(k),'.'];
    plot(om_new(k,:)/2/pi/1e3,om_new(k,:)./wavenumber,s,'LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Phase velocity');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
axis([0 5000000/1e3 0 10000]);

figure(17); hold on;
for k=1:num_of_modes
    s=[c(k),'.'];
    plot(om_new(k,:)/2/pi/1e3,cg_new(k,:),s,'LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Group velocity');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
axis([0 5000000/1e3 0 10000]);

figure(18); hold on;

plot(om_new(:,:)/2/pi/1e3,wavenumber,'LineWidth',2);

Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Wavenumber [1/m]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
axis([0 5000000/1e3 0 1800]);

%%
return
%IF3 = zeros(3,length(beta));
%mode_num = zeros(3,length(beta));
CG3 = zeros(3, length(beta));
CP3 = zeros(3, length(beta));
K3 = zeros(3, length(beta));
for j=1:length(beta)
    cc=0;
    for kk = 1:length(M)
        [aa,ii]=min(abs(FREQ(kk,:,j)-freq*1e3));
        if(ii~=1)
            cc=cc+1;
            %IF3(cc,j)=ii;
            %mode_num(cc,j)=kk;
            % interpolate profiles for selected frequency for polar plots
            CG3(cc,j) = interp1(FREQ(kk,[ii-1,ii,ii+1],j), CG(kk,[ii-1,ii,ii+1],j),freq*1e3,'linear','extrap');
            K3(cc,j) = interp1(FREQ(kk,[ii-1,ii,ii+1],j), wavenumber([ii-1,ii,ii+1]),freq*1e3,'linear','extrap');
            CP3(cc,j) = freq*1e3*2*pi/K3(cc,j);
            %CG3_(cc,j)=CG(kk,ii,j);
            if(cc==3) break; end
        end
    end
end

veloc_lim = 10000;
figure(3);
polarplot(beta*pi/180,squeeze(CG3(1,:)),'k.');
thetalim([0 90]);
hold on;
polarplot(beta*pi/180,squeeze(CG3(2,:)),'r.');
polarplot(beta*pi/180,squeeze(CG3(3,:)),'b.');
title(['Group velocity ',num2str(freq),' [kHz]']);
set(gcf,'Color','white');

figure(4);
polarplot(beta*pi/180,squeeze(CP3(1,:)),'k.');
thetalim([0 90]);
hold on;
polarplot(beta*pi/180,squeeze(CP3(2,:)),'r.');
polarplot(beta*pi/180,squeeze(CP3(3,:)),'b.');
title(['Phase velocity ',num2str(freq),' [kHz]']);
set(gcf,'Color','white');

figure(4);
polarplot(beta*pi/180,squeeze(K3(1,:)),'k.');
thetalim([0 90]);
hold on;
polarplot(beta*pi/180,squeeze(K3(2,:)),'r.');
polarplot(beta*pi/180,squeeze(K3(3,:)),'b.');
title(['Wavenumber ',num2str(freq),' [kHz]']);
set(gcf,'Color','white');

save('CFRP_4lay_0deg_1_5_mm.mat','beta','K3','CG3','CP3');