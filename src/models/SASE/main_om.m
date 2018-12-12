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
c='rgbkrbkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrgbkrbgkrbgkrgbkrbkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrgbkrbgkrbgkrgbkrbkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrbgkrgbkrbgkrbgk';
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

om_new = om_real;
ind = zeros(length(M),number_of_wavenumber_points);
cg_new=cg;
for i=1:num_of_modes
    ind(i,number_of_wavenumber_points)=i;
end
for i=1:num_of_modes
    I=i;cc=0;
for k=number_of_wavenumber_points-1:-1:4
   
    if (cc==0) % first order Taylor approximation
        cc=cc+1;
        gamma1=cg(i,k-1);
        gamma0=om_real(i,k-1);  
        om_t=gamma0-gamma1*k_step; % first order Taylor approximation (backward)
        cg_t=2*cg2-cg1; % linear Taylor approximation
    elseif (cc==1) % second order Taylor approximation
        
    end
    if(I~=i && cc==0)
        
        cc=cc+1;
        cg0=cg(i,k);
        cg1=cg(i,k-1);
        cg2=cg(I,k-2);
    
        gamma0=om_real(I,k-2);  
    elseif (I~=i && cc==1)
        cc=cc+1;
        cg0=cg(i,k);
        cg1=cg(I,k-1);
        cg2=cg(I,k-2);
    
        gamma0=om_real(I,k-2);  
    elseif (I~=i && cc>=2)
        cc=cc+1;
        cg0=cg(I,k);
        cg1=cg(I,k-1);
        cg2=cg(I,k-2);
    
        gamma0=om_real(I,k-2);
    else
        cg0=cg(i,k);
        cg1=cg(i,k-1);
        cg2=cg(i,k-2);
    
        gamma0=om_real(i,k-2);
    end
    
   
    gamma1=cg2;
    gamma2=1/2*(cg2-cg1)/k_step;
    gamma3=1/6*(cg2-2*cg1+cg0)/(k_step^2);
    
    beta2=(gamma2^2-gamma3*gamma1)/(gamma1^2-gamma2*gamma0);
    beta1=-(gamma2+gamma0*beta2)/gamma1;
    alpha0=gamma0;
    alpha1=gamma1+gamma0*beta1;
    
    %om_p=(alpha0 + alpha1 * k_step)/(1 + beta1*k_step + beta2*k_step^2); %Pade expansion of order [1/2](forward)
    om_p=(alpha0 - alpha1 * k_step)/(1 -( beta1*k_step + beta2*k_step^2)); % Pade expansion of order [1/2](backward)
    om_t=gamma0-gamma1*k_step-gamma2*k_step^2-gamma3*k_step^3; % Third order Taylor approximation (backward)
    
    cg_t=2*cg2-cg1; % linear Taylor approximation
    
    delta_om = abs((om_p-om_real(:,k-3))/om_p); % deviation of angular frequency in Pade method
    delta_cg = abs((cg_t-cg(:,k-3))/cg_t);  % deviation of group velocity 
    delta = delta_om.^2+delta_cg.^2;
    
%     delta_om = abs((om_t-om_real(:,k-3))/om_t); % deviation of angular frequency in Taylor method
%     delta_cg = abs((cg_t-cg(:,k-3))/cg_t);  % deviation of group velocity in Taylor method
%     delta = delta_om.^2+delta_cg.^2;
%     
    
    [a,I] = min(delta_om);
    om_new(i,k-3) = om_real(I,k-3);
    cg_new(i,k-3) = cg(I,k-3);
    ind(i,k-3) = I;
    
    
%     close all;
%     figure
%     plot(om_p,wavenumber(k-3),'g*');hold on;
%     plot(om_new(i,k-3),wavenumber(k-3),'bd');            
%     plot(om_real(i,k:-1:k-3),wavenumber(k:-1:k-3),'ro');
%     figure  
%     plot(cg_t,wavenumber(k-3),'g*');hold on;
%     plot(cg_new(i,k-3),wavenumber(k-3),'bd');            
%     plot(cg(i,k:-1:k-3),wavenumber(k:-1:k-3),'ro');
%     
%     pause;
end
end
%%
figure(16); hold on;
for k=1:num_of_modes
    s=[c(k),'.'];
    plot(om_new(k,:)/2/pi/1e3,om_new(k,:)./wavenumber,s,'LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(beta(1)),' [deg]']);
axis([0 5000000/1e3 0 10000]);


return;
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
        [row,col]=size(wave_number_imag);
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

figure(13); hold on;
for k=1:num_of_modes
    plot(f,wave_number_real_sort(k,:),'k.','LineWidth',2);
end
Hxl=xlabel('Frequency [kHz]');
Hyl=ylabel('Wavenumber [1/m]');
set(Hxl,'FontSize',12);set(Hyl,'FontSize',12);
set(gcf,'Color','white');set(gca,'FontSize',11);
title([num2str(theta),' [deg]']);

save([fname,'.mat'],'f','wave_number_real_sort');