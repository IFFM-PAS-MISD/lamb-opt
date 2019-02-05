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
beta = 0:15:90; % angles for dispersion curves in polar plot for fixed frequency
np = 3;
nele_layer = 1;
%freq=200; % frequency for polar plots [kHz]
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
    %wavenumber(k) = 100;
    
    [cg(:,k),mode_shapes, om_real(:,k), om_imag(:,k)] = SASE_om(K,M,beta(1),wavenumber(k));
    %[cg2(:,k),mode_shapes, om_real2(:,k), om_imag2(:,k)] = SASE_om(K,M,beta(1),wavenumber(k)+wavenumber_min);
end


%%
c='bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';
figure(12); hold on;
for k=1:num_of_modes
    s=[c(k),'.'];
    plot(om_real(k,:)/2/pi/1e3,wavenumber,s,'LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
%Hxl=xlabel('Frequency [Hz]');
Hyl=ylabel('Wavenumber [1/m]');

set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
axis([0 500000/1e3 0 1800]);


figure(13); hold on;
for k=1:num_of_modes
    s=[c(k),'.'];
    plot(om_real(k,:)/2/pi/1e3,om_real(k,:)./wavenumber,s,'LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Phase velocity');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
axis([0 500000/1e3 0 10000]);

figure(14); hold on;
for k=1:num_of_modes
    s=[c(k),'.'];
    plot(om_real(k,:)/2/pi/1e3,(real(cg(k,:))),s,'LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Group velocity');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
axis([0 500000/1e3 0 10000]);

% figure(15); hold on;
% for k=1:num_of_modes
%     s=[c(k),'.'];
%     plot(om_real(k,:)/2/pi/1e3,(om_real2(k,:)-om_real(k,:))./wavenumber_min,s,'LineWidth',2);
% end
% Hxl=xlabel('Frequency [kHz]');
% Hyl=ylabel('Group velocity');
% set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
% set(gcf,'Color','white');set(gca,'FontSize',11);
% title([num2str(beta(1)),' [deg]']);
% axis([0 500000/1e3 0 10000]);

%% mode-tracing
% Pade expansion method
[cg_new,om_new] = mode_tracing(cg,om_real,k_step);

%%

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
return;

save([fname,'.mat'],'f','wave_number_real_sort');