% genetic algorithm - identification of elastic constants in fibre reinforced composite laminate

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true
number_of_modes_considered = 4; % number of modes considered in calculation of objective function score
%% Load parameters which are used in experiment
% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );
% filename of parameter data
 filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_param'; 
 load([data_path,filename]); % wavenumber_max fmax beta number_of_wavenumber_points
% load experimental data file
exp_filename = {'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF'};
load([data_path,exp_filename{1}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points 
%% Input for SASE
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
layup = [0 0 0 0 0 0 0 0];
nlayers = length(layup);
h = [zeros(nlayers,1)+1]* 3e-3/nlayers; % thickness of layers;
% Stacking direction
stack_dir = 1;

 %% Plot best case
radians = false;
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 8; 
[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);
load project_paths projectroot src_path;
run([src_path,filesep,'models',filesep,'SASE',filesep,'inputs',filesep,'Fabric_1.m']);
%%
load('test7_1');
MAXGEN = 40;
rho_m1 = PBest(MAXGEN,1);
rho_f1 = PBest(MAXGEN,2);
e11_m1 = PBest(MAXGEN,3)/1e9;
e11_f1 = PBest(MAXGEN,4)/1e9;
ni12_m1 = PBest(MAXGEN,5);
ni12_f1 = PBest(MAXGEN,6);
vol_01 = PBest(MAXGEN,7);
e22_f1 = 0.1*e11_f;
ni23_f1 = PBest(MAXGEN,6);
%% Mechanical properties  
 [C11_1,C12_1,C13_1,C21_1,C22_1,C23_1,C31_1,C32_1,C33_1,C44_1,C55_1,C66_1,rho_1] = ...
        compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_01, ...
        e11_m1, ni12_m1, rho_m1, e11_f1, e22_f1, ni12_f1, ni23_f1, rho_f1,false);
%% SASE
[wavenumber,CG,FREQ_1] = main_SASE(rho_1,C11_1,C12_1,C13_1,C22_1,C23_1,C33_1,C44_1,C55_1,C66_1,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
%%
load('test7_2');
rho_m2 = PBest(MAXGEN,1);
rho_f2 = PBest(MAXGEN,2);
e11_m2 = PBest(MAXGEN,3)/1e9;
e11_f2 = PBest(MAXGEN,4)/1e9;
ni12_m2 = PBest(MAXGEN,5);
ni12_f2 = PBest(MAXGEN,6);
vol_02 = PBest(MAXGEN,7);
e22_f2 = 0.1*e11_f;
ni23_f2 = PBest(MAXGEN,6);

%% Mechanical properties  
 [C11_2,C12_2,C13_2,C21_2,C22_2,C23_2,C31_2,C32_2,C33_2,C44_2,C55_2,C66_2,rho_2] = ...
        compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_02, ...
        e11_m2, ni12_m2, rho_m2, e11_f2, e22_f2, ni12_f2, ni23_f2, rho_f2,false);
%% SASE
[wavenumber,CG,FREQ_2] = main_SASE(rho_2,C11_2,C12_2,C13_2,C22_2,C23_2,C33_2,C44_2,C55_2,C66_2,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);

for j=1:number_of_angles
    figure(j)
        set(gcf,'Color','w');
        %imagesc(fvec(2:end)/1e3, wavenumber(2:end,j)/2, squeeze(abs(Data_polar(j,2:end,2:end)))); 
        imagesc(fvec(2:end)/1e3, wavenumber(2:end,j), squeeze(abs(Data_polar(j,2:end,2:end)))); 
        set(gca,'YDir','normal'); 
        %axis([0 600 0 2000]);
        %%axis([0 350 0 2000]);
        axis([0 600 0 min(wavenumber_max)]);
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
        if(radians)
            ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
        else
            ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
        end
        colormap jet; 
        %axis tight; 
        caxis([0 max(caxis)/3]); 

        hold on;
        fvec1_1=squeeze(FREQ_1(1,:,j)); % mode 1, angle j
        fvec2_1=squeeze(FREQ_1(2,:,j)); % mode 2, angle j
        fvec3_1=squeeze(FREQ_1(3,:,j)); % mode 3, angle j
        fvec4_1=squeeze(FREQ_1(4,:,j)); % mode 4, angle j
        fvec5_1=squeeze(FREQ_1(5,:,j)); % mode 1, angle j
        fvec6_1=squeeze(FREQ_1(6,:,j)); % mode 2, angle j
        fvec7_1=squeeze(FREQ_1(7,:,j)); % mode 3, angle j
        fvec8_1=squeeze(FREQ_1(8,:,j)); % mode 4, angle j
        kvec=squeeze(wavenumber(:,j)); % angle j
        plot(fvec1_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec2_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec3_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec4_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec5_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec6_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec7_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec8_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        
        fvec1_2=squeeze(FREQ_2(1,:,j)); % mode 1, angle j
        fvec2_2=squeeze(FREQ_2(2,:,j)); % mode 2, angle j
        fvec3_2=squeeze(FREQ_2(3,:,j)); % mode 3, angle j
        fvec4_2=squeeze(FREQ_2(4,:,j)); % mode 4, angle j
        fvec5_2=squeeze(FREQ_2(5,:,j)); % mode 1, angle j
        fvec6_2=squeeze(FREQ_2(6,:,j)); % mode 2, angle j
        fvec7_2=squeeze(FREQ_2(7,:,j)); % mode 3, angle j
        fvec8_2=squeeze(FREQ_2(8,:,j)); % mode 4, angle j
        kvec=squeeze(wavenumber(:,j)); % angle j
        plot(fvec1_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec2_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec3_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec4_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec5_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec6_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec7_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec8_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        % get(gca,'FontName'); % default 'Helvetica'
        %set(gca,'FontName','Arial');
        %set(gca,'FontName','Helvetica');
        set(gca,'FontName','Times');
        fig = gcf;
        %title(['angle ', num2str(beta(j)),' deg'],'Fontsize',12);
        %title(['angle ', num2str(beta(j))],'Fontsize',12);
        %title(['angle ', num2str(beta(j)),' deg'],'Fontsize',12,'interpreter','latex');
        title({[num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
end  

for j=1:number_of_angles
    figure(number_of_angles+j)
        set(gcf,'Color','w');
        
        set(gca,'YDir','normal'); 
       
        axis([0 600 0 min(wavenumber_max)]);
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
        if(radians)
            ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
        else
            ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
        end
        colormap jet; 
        %axis tight; 
        caxis([0 max(caxis)/3]); 

        hold on;
        fvec1_1=squeeze(FREQ_1(1,:,j)); % mode 1, angle j
        fvec2_1=squeeze(FREQ_1(2,:,j)); % mode 2, angle j
        fvec3_1=squeeze(FREQ_1(3,:,j)); % mode 3, angle j
        fvec4_1=squeeze(FREQ_1(4,:,j)); % mode 4, angle j
        fvec5_1=squeeze(FREQ_1(5,:,j)); % mode 1, angle j
        fvec6_1=squeeze(FREQ_1(6,:,j)); % mode 2, angle j
        fvec7_1=squeeze(FREQ_1(7,:,j)); % mode 3, angle j
        fvec8_1=squeeze(FREQ_1(8,:,j)); % mode 4, angle j
        kvec=squeeze(wavenumber(:,j)); % angle j
        plot(fvec1_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec2_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec3_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec4_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec5_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec6_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec7_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec8_1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        
        fvec1_2=squeeze(FREQ_2(1,:,j)); % mode 1, angle j
        fvec2_2=squeeze(FREQ_2(2,:,j)); % mode 2, angle j
        fvec3_2=squeeze(FREQ_2(3,:,j)); % mode 3, angle j
        fvec4_2=squeeze(FREQ_2(4,:,j)); % mode 4, angle j
        fvec5_2=squeeze(FREQ_2(5,:,j)); % mode 1, angle j
        fvec6_2=squeeze(FREQ_2(6,:,j)); % mode 2, angle j
        fvec7_2=squeeze(FREQ_2(7,:,j)); % mode 3, angle j
        fvec8_2=squeeze(FREQ_2(8,:,j)); % mode 4, angle j
        kvec=squeeze(wavenumber(:,j)); % angle j
        plot(fvec1_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec2_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec3_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec4_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec5_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec6_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec7_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        plot(fvec8_2(2:end)/1e3,kvec(2:end),'m','linewidth',1);
        % get(gca,'FontName'); % default 'Helvetica'
        %set(gca,'FontName','Arial');
        %set(gca,'FontName','Helvetica');
        set(gca,'FontName','Times');
        fig = gcf;
        %title(['angle ', num2str(beta(j)),' deg'],'Fontsize',12);
        %title(['angle ', num2str(beta(j))],'Fontsize',12);
        %title(['angle ', num2str(beta(j)),' deg'],'Fontsize',12,'interpreter','latex');
        title({[num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
end 