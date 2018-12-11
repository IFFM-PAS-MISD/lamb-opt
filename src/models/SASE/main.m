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
beta = 0:5:90;
np = 3;
nele_layer = 1;
freq=200; % frequency for polar plots [kHz]
theta=0; % angle for dispersion curves
fmin=1; %minimal frequency [kHz]
frmax=20; % maximal frequency for dispersion curves [kHz]
number_of_freq_points=500;
%number_of_freq_points=750;
%number_of_freq_points=1000;
fr_step=(frmax-fmin)/(number_of_freq_points-1); % frequency step [kHz]
nbeta = length(beta);
delta_freq=fr_step*1/100; % fr_step*10 [Hz]
num_of_modes=9; % number of modes to display on dispersion curves
type='composite'; % type='composite' or type='isotropic'
%type='isotropic';
%fname='Tomek_carbon_epoxy_200kHz_s'; % filename for saving figures
%fname='alu_1mm_25C'; % filename for saving figures
%fname='alu_3mm_25C'; % filename for saving figures
%fname='my_45mm_fibers'; % filename for saving figures
%fname='my_45mm_fibers_vol20'; % filename for saving figures
%fname='my_45mm_fibers_vol10'
%fname='PANmodel'; % filename for saving figures
fname='Airbus';
sorting = 'no';
%% SAFE
disp('SASE...');
freq2=freq+delta_freq;
om=2*pi*freq;
om2=2*pi*freq2;
[mode_shapes, wave_number_real, wave_number_imag] = SASE(nele_layer,np,beta, freq);
[mode_shapes2, wave_number_real2, wave_number_imag2] = SASE(nele_layer,np,beta, freq2);
%% Roots sorting
switch sorting
    case 'yes'
[row,col]=size(wave_number_imag);
wave_number_imag_sort=zeros(row,col);
wave_number_real_sort=zeros(row,col);
wave_number_imag_sort2=zeros(row,col);
wave_number_real_sort2=zeros(row,col);
switch type
    case 'composite'
        [Ind]=roots_sorting(wave_number_imag,beta);
        [Ind2]=roots_sorting(wave_number_imag2,beta);
    case 'isotropic'
        [Ind]=roots_sorting(wave_number_real,beta);
        [Ind2]=roots_sorting(wave_number_real2,beta);
end
for i=1:row
    for j=1:col
        wave_number_imag_sort(i,j) = wave_number_imag(Ind(i,j),j);
        wave_number_real_sort(i,j) = wave_number_real(Ind(i,j),j);
        wave_number_imag_sort2(i,j) = wave_number_imag2(Ind2(i,j),j);
        wave_number_real_sort2(i,j) = wave_number_real2(Ind2(i,j),j);
    end
end
    case 'no'
        wave_number_imag_sort = wave_number_imag;
        wave_number_real_sort =  wave_number_real;
        wave_number_imag_sort2 = wave_number_imag2;
        wave_number_real_sort2 =  wave_number_real2;
end
%% Polar plot

% in Np/m -- convert to dB/m -- 1 Np/m = 20/log(10) dB/m

attenuation_mod1 = wave_number_imag_sort(1,:);
% attenuation_mod1 = attenuation_mod1 * 20/log(10);
cp1=om./wave_number_real_sort(1,:);
%cp2=om2./wave_number_real_sort2(1,:);
phase_veloc_mod1=cp1; %[km/s]
group_veloc_mod1=(om2-om)./(wave_number_real_sort2(1,:)-wave_number_real_sort(1,:)); %[km/s]
%group_veloc_mod1=cp2.^2./(cp2-(om2.*(cp2-cp1)./(om2-om))); %[km/s]

attenuation_mod2 = wave_number_imag_sort(2,:);
% attenuation_mod2 = attenuation_mod2 * 20/log(10);
cp1=om./wave_number_real_sort(2,:);
%cp2=om2./wave_number_real_sort2(2,:);
phase_veloc_mod2=cp1; %[km/s]
group_veloc_mod2=(om2-om)./(wave_number_real_sort2(2,:)-wave_number_real_sort(2,:)); %[km/s]
%group_veloc_mod2=cp2.^2./(cp2-(om2.*(cp2-cp1)./(om2-om))); %[km/s]

attenuation_mod3 = wave_number_imag_sort(3,:);
% attenuation_mod3 = attenuation_mod3 * 20/log(10);
cp1=om./wave_number_real_sort(3,:);
%cp2=om2./wave_number_real_sort2(3,:);
phase_veloc_mod3=cp1; %[km/s]

group_veloc_mod3=(om2-om)./(wave_number_real_sort2(3,:)-wave_number_real_sort(3,:)); %[km/s]
%group_veloc_mod3=cp2.^2./(cp2-(om2.*(cp2-cp1)./(om2-om))); %[km/s]


att_lim = 1.0*max([max(attenuation_mod1),max(attenuation_mod2),max(attenuation_mod3)]);
phase_veloc_lim = 1.0*max([max(phase_veloc_mod1),max(phase_veloc_mod2),max(phase_veloc_mod3)]);
group_veloc_lim = 1.0*max([max(group_veloc_mod1),max(group_veloc_mod2),max(group_veloc_mod3)]);
%%
% attenuation plots

figure(1);
newPolar(beta*pi/180,attenuation_mod1,att_lim,'k.');
hold on;
newPolar(beta*pi/180,attenuation_mod2,att_lim,'r.');
newPolar(beta*pi/180,attenuation_mod3,att_lim,'b.');
title(['Attenuation ',num2str(freq),' [kHz]']);
set(gcf,'Color','white');
%print(gcf,'-dtiffn','-r600',['figs\Attenuation_polar_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Attenuation_polar_',fname,'.png']);
saveas(gcf,['figs\Attenuation_polar_',fname,'.fig'], 'fig');

figure(2); hold on;
% plot(beta,attenuation_mod1,'k.','LineWidth',2);
% plot(beta,attenuation_mod2,'r.','LineWidth',2);
% plot(beta,attenuation_mod3,'b.','LineWidth',2);
plot(beta,attenuation_mod1,'k.','LineWidth',2);
plot(beta,attenuation_mod2,'k.','LineWidth',2);
plot(beta,attenuation_mod3,'k.','LineWidth',2);

axis([0 360 0 att_lim])
%legend('mod1','mod2','mod3');
Hxl=xlabel('Angle');
Hyl=ylabel('Attenuation [Np/m]');
title([num2str(freq),' [kHz]']);
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
%print(gcf,'-dtiffn','-r600',['figs\Attenuation_angle_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Attenuation_angle_',fname,'.png']);
saveas(gcf,['figs\Attenuation_angle_',fname,'.fig'], 'fig');
% phase velocities plots
figure(3);
newPolar(beta*pi/180,phase_veloc_mod1,phase_veloc_lim,'k.');
hold on;
newPolar(beta*pi/180,phase_veloc_mod2,phase_veloc_lim,'r.');
newPolar(beta*pi/180,phase_veloc_mod3,phase_veloc_lim,'b.');
title(['Phase velocity ',num2str(freq),' [kHz]']);
set(gcf,'Color','white');
%legend('mod2','mod1','mod3');
%print(gcf,'-dtiffn','-r600',['figs\Phase_velocity_polar_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Phase_velocity_polar_',fname,'.png']);
saveas(gcf,['figs\Phase_velocity_polar_',fname,'.fig'], 'fig');

figure(4); hold on;
plot(beta,phase_veloc_mod1,'k.','LineWidth',2);
plot(beta,phase_veloc_mod2,'r.','LineWidth',2);
plot(beta,phase_veloc_mod3,'b.','LineWidth',2);
axis([0 360 0 phase_veloc_lim])
%legend('mod1','mod2','mod3');
%legend('mod2','mod1','mod3');
Hxl=xlabel('Angle');
Hyl=ylabel('Phase velocity');
title([num2str(freq),' [kHz]']);
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
%print(gcf,'-dtiffn','-r600',['figs\Phase_velocity_angle_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Phase_velocity_angle_',fname,'.png']);
saveas(gcf,['figs\Phase_velocity_angle_',fname,'.fig'], 'fig');

% group velocities plots
figure(5);
% newPolar(beta*pi/180,group_veloc_mod1,group_veloc_lim,'k.');
% hold on;
% newPolar(beta*pi/180,group_veloc_mod2,group_veloc_lim,'r.');
% newPolar(beta*pi/180,group_veloc_mod3,group_veloc_lim,'b.');
newPolar(beta*pi/180,group_veloc_mod1,group_veloc_lim,'k.');
hold on;
newPolar(beta*pi/180,group_veloc_mod2,group_veloc_lim,'k.');
newPolar(beta*pi/180,group_veloc_mod3,group_veloc_lim,'k.');
title(['Group velocity ',num2str(freq),' [kHz]']);
set(gcf,'Color','white');
%print(gcf,'-dtiffn','-r600',['figs\Group_velocity_polar_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Group_velocity_polar_',fname,'.png']);
saveas(gcf,['figs\Group_velocity_polar_',fname,'.fig'], 'fig');

figure(6); hold on;
plot(beta,group_veloc_mod1,'k.','LineWidth',2);
plot(beta,group_veloc_mod2,'r.','LineWidth',2);
plot(beta,group_veloc_mod3,'b.','LineWidth',2);
axis([0 360 0 group_veloc_lim]);
%legend('mod1','mod2','mod3');
Hxl=xlabel('Angle');
Hyl=ylabel('Group velocity');
title([num2str(freq),' [kHz]']);
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
%print(gcf,'-dtiffn','-r600',['figs\Group_velocity_angle_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Group_velocity_angle_',fname,'.png']);
saveas(gcf,['figs\Group_velocity_angle_',fname,'.fig'], 'fig');
%return;
%% dispersion curves
frmin=fmin; %
frn=round((frmax-frmin)/fr_step);
wave_number_real = zeros(row,frn);
wave_number_real2 = zeros(row,frn);
wave_number_imag = zeros(row,frn);
wave_number_imag2 = zeros(row,frn);


f=zeros(1,frn);
om=zeros(1,frn);
% in Np/m -- convert to dB/m -- 1 Np/m = 20/log(10) dB/m
disp('SASE...');
for k=1:frn
    [k frn]
    freq=frmin+(k-1)*fr_step;
    f(1,k)=freq;
    %% SAFE
    om(k)=2*pi*freq;
    freq2=freq+delta_freq;
    om2(k)=2*pi*freq2;
    
    [mode_shapes(:,k), wave_number_real(:,k), wave_number_imag(:,k)] = SASE(nele_layer,np,theta, freq);
    [mode_shapes2(:,k), wave_number_real2(:,k), wave_number_imag2(:,k)] = SASE(nele_layer,np,theta, freq2);
end

%% Roots sorting
switch sorting
    case 'yes'
    [row,col]=size(wave_number_imag);
    wave_number_imag_sort=zeros(row,frn);
    wave_number_real_sort=zeros(row,frn);
    wave_number_imag_sort2=zeros(row,frn);
    wave_number_real_sort2=zeros(row,frn);

switch type
    case 'composite'
        [Ind]=roots_sorting(wave_number_imag,f);
        [Ind2]=roots_sorting(wave_number_imag2,f);
    case 'isotropic'
        [Ind]=roots_sorting(wave_number_real,f);
        [Ind2]=roots_sorting(wave_number_real2,f);
end
    for i=1:row
        for j=1:col
            wave_number_imag_sort(i,j) = wave_number_imag(Ind(i,j),j);
            wave_number_real_sort(i,j) = wave_number_real(Ind(i,j),j);
            wave_number_imag_sort2(i,j) = wave_number_imag2(Ind2(i,j),j);
            wave_number_real_sort2(i,j) = wave_number_real2(Ind2(i,j),j);
        end
    end
    case 'no'
         wave_number_imag_sort = wave_number_imag;
         wave_number_real_sort =  wave_number_real;
         wave_number_imag_sort2 = wave_number_imag2;
         wave_number_real_sort2 =  wave_number_real2;
end
   %%

    attenuation_mod1_fr = wave_number_imag_sort(1,:);
    % attenuation_mod1 = attenuation_mod1 * 20/log(10);
    phase_veloc_mod1_fr=om./wave_number_real_sort(1,:); %[km/s]
    group_veloc_mod1_fr=(om2-om)./(wave_number_real_sort2(1,:)-wave_number_real_sort(1,:)); %[km/s]
    %group_veloc_mod1_fr=phase_veloc_mod1_fr(2:end).^2./(phase_veloc_mod1_fr(2:end)-om(2:end).*(phase_veloc_mod1_fr(2:end)-phase_veloc_mod1_fr(1:end-1))./(om(2:end)-om(1:end-1))); %[km/s]

    attenuation_mod2_fr = wave_number_imag_sort(2,:);
    % attenuation_mod2 = attenuation_mod2 * 20/log(10);
    phase_veloc_mod2_fr=om./wave_number_real_sort(2,:); %[km/s]
    group_veloc_mod2_fr=(om2-om)./(wave_number_real_sort2(2,:)-wave_number_real_sort(2,:)); %[km/s]
    %group_veloc_mod2_fr=phase_veloc_mod2_fr(2:end).^2./(phase_veloc_mod2_fr(2:end)-om(2:end).*(phase_veloc_mod2_fr(2:end)-phase_veloc_mod2_fr(1:end-1))./(om(2:end)-om(1:end-1))); %[km/s]

    attenuation_mod3_fr = wave_number_imag_sort(3,:);
    % attenuation_mod3 = attenuation_mod3 * 20/log(10);
    phase_veloc_mod3_fr=om./wave_number_real_sort(3,:); %[km/s]
    group_veloc_mod3_fr=(om2-om)./(wave_number_real_sort2(3,:)-wave_number_real_sort(3,:)); %[km/s]
    %group_veloc_mod3_fr=phase_veloc_mod3_fr(2:end).^2./(phase_veloc_mod3_fr(2:end)-om(2:end).*(phase_veloc_mod3_fr(2:end)-phase_veloc_mod3_fr(1:end-1))./(om(2:end)-om(1:end-1))); %[km/s]

% plots

figure(7); hold on;
plot(f,attenuation_mod1_fr,'k.','LineWidth',2);
plot(f,attenuation_mod2_fr,'r.','LineWidth',2);
plot(f,attenuation_mod3_fr,'b.','LineWidth',2);
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Attenuation [Np/m]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(theta),' [deg]']);
%print(gcf,'-dtiffn','-r600',['figs\Attenuation_dispersion_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Attenuation_dispersion_',fname,'.png']);
saveas(gcf,['figs\Attenuation_dispersion_',fname,'.fig'], 'fig');

figure(8); hold on;
plot(f,phase_veloc_mod1_fr,'k.','LineWidth',2);
plot(f,phase_veloc_mod2_fr,'r.','LineWidth',2);
plot(f,phase_veloc_mod3_fr,'b.','LineWidth',2);
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Phase velocity [km/s]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
% axis([0 frmax 0 10]);
axis([0 frmax 0 5]);
%title([num2str(theta),' [deg]']);
%print(gcf,'-dtiffn','-r600',['figs\Phase_velocity_dispersion_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Phase_velocity_dispersion_',fname,'.png']);
saveas(gcf,['figs\Phase_velocity_dispersion_',fname,'.fig'], 'fig');
phase_disp=[f', phase_veloc_mod1_fr'*1e3, phase_veloc_mod2_fr'*1e3,phase_veloc_mod3_fr'*1e3];
save(['phase_disp_',fname,'.txt'],'phase_disp','-ascii');

figure(9); hold on;
plot(f,group_veloc_mod1_fr,'k.','LineWidth',2);
plot(f,group_veloc_mod2_fr,'r.','LineWidth',2);
plot(f,group_veloc_mod3_fr,'b.','LineWidth',2);
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Group velocity [km/s]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
% axis([0 frmax 0 10]);
axis([0 frmax 0 5]);
%title([num2str(theta),' [deg]']);
%print(gcf,'-dtiffn','-r600',['figs\Group_velocity_dispersion_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Group_velocity_dispersion_',fname,'.png']);
saveas(gcf,['figs\Group_velocity_dispersion_',fname,'.fig'], 'fig');
group_disp=[f', group_veloc_mod1_fr'*1e3, group_veloc_mod2_fr'*1e3,group_veloc_mod3_fr'*1e3];
save(['group_disp_',fname,'.txt'],'group_disp','-ascii');

figure(10); hold on;
for k=1:num_of_modes
plot(f,wave_number_imag_sort(k,:),'k.','LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Attenuation [Np/m]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(theta),' [deg]']);
%print(gcf,'-dtiffn','-r600',['figs\Atteniuation_dispersion_all_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Atteniuation_dispersion_all_',fname,'.png']);
saveas(gcf,['figs\Atteniuation_dispersion_all_',fname,'.fig'], 'fig');

figure(11); hold on;
for k=1:num_of_modes
    cp=om./wave_number_real_sort(k,:);
    plot(f,cp,'k.','LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Phase velocity [km/s]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
axis([0 frmax 0 40]);
title([num2str(theta),' [deg]']);
%print(gcf,'-dtiffn','-r600',['figs\Phase_velocity_dispersion_all_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Phase_velocity_dispersion_all_',fname,'.png']);
saveas(gcf,['figs\Phase_velocity_dispersion_all_',fname,'.fig'], 'fig');

figure(12); hold on;
for k=1:num_of_modes
     cp=om./wave_number_real_sort(k,:);
     cg=cp(1:end-1).^2./(cp(1:end-1)-om(1:end-1).*(cp(2:end)-cp(1:end-1))./(om(2:end)-om(1:end-1))); %[km/s]    
     plot(f(1:end-1),cg,'k.','LineWidth',2);   
 %   cg=(om2-om)./(wave_number_real_sort2(k,:)-wave_number_real_sort(k,:)); %[km/s]
  %  plot(f,cg,'k.','LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Group velocity [km/s]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
axis([0 frmax 0 10]);
title([num2str(theta),' [deg]']);
%print(gcf,'-dtiffn','-r600',['figs\Group_velocity_dispersion_all_',fname,'.tiff']);
print(gcf,'-dpng','-r300',['figs\Group_velocity_dispersion_all_',fname,'.png']);
saveas(gcf,['figs\Group_velocity_dispersion_all_',fname,'.fig'], 'fig');

save([fname,'.mat'],'f','wave_number_real_sort');