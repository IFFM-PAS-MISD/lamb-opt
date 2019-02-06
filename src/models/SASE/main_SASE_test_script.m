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
beta = 0:15:90; % angles for dispersion curves in polar plot
np = 3;
nele_layer = 1;

wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]

wavenumber_max = [3657.8,3786.8,4223.7,5172.9,4223.7,3786.8,3657.8]'; % maximal wavenumber for dispersion curves [1/m]
number_of_wavenumber_points=512;


%% SAFE
disp('SASE...');

%
%% Input

rhom = 1250; % kg/m^3
rhof = 1900; % kg/m^3
em = 3.43e9; % Pa
ef = 240e9; % Pa
nim = 0.35;
nif =  0.2; 
vol = 0.5;
layup = [0 90 90 0];
nlayers = length(layup);
h = [1,1,1,1]'* 1.5e-3/nlayers; % thickness of layers;
% Stacking direction
stack_dir = 1;
   

% homogenization of unidirectional fibre reinforce composite
[rho,e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23] = homogenization(rhom,rhof,em,ef,nim,nif,vol);
% Elastic constants of composite lamina in terms of principal material directions
[C11,C12,C13,C22,C23,C33,C44,C55,C66] = lamina_3D(e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23);

[wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);

% identify A0 and S0 mode
FREQ_A0 = zeros(number_of_wavenumber_points,length(beta));
FREQ_S0 = zeros(number_of_wavenumber_points,length(beta));
CG_A0 = zeros(number_of_wavenumber_points,length(beta));
CG_S0 = zeros(number_of_wavenumber_points,length(beta));
[num_of_modes,number_of_wavenumber_points,number_of_angles] = size(CG);
for j=1:number_of_angles
    f1=FREQ(:,2,j);
    I3=find(f1<50e3); % should find 3 modes existing in range 0-50 kHz, namely A0, S0 and SH0
    [~,A0_ind] = min(CG(I3,4,j)); % A0 mode index
    [~,S0_ind] = max(CG(I3,4,j)); % S0 mode index
    FREQ_A0(:,j)=FREQ(I3(A0_ind),:,j);
    CG_A0(:,j)=CG(I3(A0_ind),:,j);
    FREQ_S0(:,j)=FREQ(I3(S0_ind),:,j);
    CG_S0(:,j)=CG(I3(S0_ind),:,j);
end
% sanity check
figure(1)
hold on;
for j=1:number_of_angles
    plot(FREQ_A0(2:end,j)/1e3,CG_A0(2:end,j),'LineWidth',2);% A0
    %plot(FREQ_S0(:,j)/1e3,CG_S0(:,j),'LineWidth',2);% S0
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Group velocity [m/s]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
figure(2)
hold on;
for j=1:number_of_angles
    %plot(FREQ_A0(:,j)/1e3,CG_A0(:,j),'LineWidth',2);% A0
    plot(FREQ_S0(2:end,j)/1e3,CG_S0(2:end,j),'LineWidth',2);% S0
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Group velocity [m/s]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
figure(3)
hold on;
for j=1:number_of_angles
    plot(FREQ_A0(2:end,j)/1e3,wavenumber(2:end,j),'LineWidth',2);% A0
    plot(FREQ_S0(2:end,j)/1e3,wavenumber(2:end,j),'LineWidth',2);% S0
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Wavenumber [1/m]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12)

return;
figure(5); hold on;

%plot(FREQ(:,:,1)/1e3,wavenumber,'LineWidth',2);
plot(FREQ(1,:,1)/1e3,wavenumber,'LineWidth',2);% A0
hold on;
plot(FREQ(4,:,1)/1e3,wavenumber,'LineWidth',2); % S0

Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Wavenumber [1/m]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
%axis([0 5000000/1e3 0 1800]);
%axis([0 500000/1e3 0 1800]);

figure(6); 
plot(FREQ(1,:,1)/1e3,CG(1,:,1),'LineWidth',2);% A0
hold on;
plot(FREQ(2,:,1)/1e3,CG(2,:,1),'LineWidth',2);

plot(FREQ(3,:,1)/1e3,CG(3,:,1),'LineWidth',2);
plot(FREQ(4,:,1)/1e3,CG(4,:,1),'LineWidth',2); % S0 0 deg
plot(FREQ(4,:,1)/1e3,CG(4,:,3),'LineWidth',2); % S0 30 deg
plot(FREQ(4,:,7)/1e3,CG(4,:,7),'LineWidth',2); % S0 90 deg


Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Group velocity [m/s]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
%axis([0 500000/1e3 0 10000]);
return;
