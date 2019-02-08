%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    main                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calls main_SASE function for a given wavenumbers (in 1/m) 
% and gets the real frequencies and group velocities.
% Dispersion curves are calculated at different angles (beta) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;
% clc;

%% Input for SASE

% Mode, frequency and angles and order of sprectal elements
beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
np = 4; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1)+1; % minimal wavenumbers [1/m]
wavenumber_max = [3657.8,3786.8,4223.7,5172.9,4223.7,3786.8,3657.8]'; % maximal wavenumber for dispersion curves [1/m]
number_of_wavenumber_points=512;
%number_of_wavenumber_points=1024;
%
%% Input for material properties

rhom = 1250; % kg/m^3
rhof = 1900; % kg/m^3
em = 3.43e9; % Pa
ef = 0.8*240e9; % Pa
nim = 0.35;
nif =  0.2; 
vol = 0.8*0.5;
layup = [0 90 0 90 90 0 90 0];
%layup = [0 90 0 90 0 90 0 90];
nlayers = length(layup);
h = [1,1,1,1]'* 3e-3/nlayers; % thickness of layers;
% Stacking direction
stack_dir = 1;
%% Mechanical properties  
% homogenization of unidirectional fibre reinforce composite
[rho,e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23] = homogenization(rhom,rhof,em,ef,nim,nif,vol);
% Elastic constants of composite lamina in terms of principal material directions
[C11,C12,C13,C22,C23,C33,C44,C55,C66] = lamina_3D(e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23);
%% SASE
disp('SASE...');
[wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
[number_of_modes,number_of_wavenumber_points,number_of_angles] = size(CG);
% for k=1:number_of_modes
% CG(k,2:end,:) = squeeze(FREQ(k,2:end,:) - FREQ(k,1:end-1,:))./(wavenumber(2:end,:)-wavenumber(1:end-1,:));
% end
% identify A0 and S0 mode
[FREQ_A0,FREQ_S0,CG_A0,CG_S0] = identify_A0_S0_modes2(FREQ,CG);

% sanity check
% all angles - fundamental modes
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
    plot(FREQ_S0(2:end,j)/1e3,real(CG_S0(2:end,j)),'LineWidth',2);% S0
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


% all modes - 0 deg
jj=1;
figure(4); hold on;
for k=1:number_of_modes
    plot(FREQ(k,2:end,jj)/1e3,wavenumber(2:end,jj),'LineWidth',2);% 
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Wavenumber [1/m]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
%axis([0 5000000/1e3 0 1800]);
%axis([0 500000/1e3 0 1800]);

figure(5); hold on;
for k=1:number_of_modes
    plot(FREQ(k,2:end,jj)/1e3,CG(k,2:end,jj),'LineWidth',2);% 
end

Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Group velocity [m/s]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
%axis([0 500000/1e3 0 10000]);

