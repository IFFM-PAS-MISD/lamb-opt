% plot numerical dispersion curves

clear all; close all;

%set(0,'defaulttextinterpreter','none');

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true
number_of_modes_considered = 4; % number of modes considered in calculation of objective function score
% figure parameters
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
%fig_width = 7; fig_height = 5; 
% create path to the numerical model data folder
modelfolder = 'genetic_algorithm';
modelname = 'ga_plain_weave_known_mass';
radians = false;
% create output path
output_path = prepare_figure_paths(modelfolder,modelname);

% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname,'_out'], filesep );
% create path to the experimental processed data folder
exp_input_data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
% and conversion to polar coordinate system
input_file = 2; % experimental data
% filename of parameter data
 filename = {'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_param',...
                    'polar_interim_499x499p_chp200_x40_18Vpp_250Hz_KXKYF_param'}; 
 load([exp_input_data_path,filename{input_file}]); % wavenumber_max fmax beta number_of_wavenumber_points
% load experimental data file
exp_filename = {'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF',... % 1 old plain weave
                            'polar_interim_499x499p_chp200_x40_18Vpp_250Hz_KXKYF'};         % 2 new plain weave
load([exp_input_data_path,exp_filename{input_file}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points  

%% Input for SASE
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
%ht = 3/1000; % [m] laminate total thickness; old plain weave
ht = 3.9/1000; % [m] laminate total thickness; new plain weave
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
layup = [0 0 0 0 0 0 0 0];
%layup =[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
nlayers = length(layup);

h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers; new plain weave
% Stacking direction
stack_dir = 1;
%%
[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);

% if(~radians)
%    wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
% end
fprintf('Making figures: %s\n', modelname);
%%
%% input material properties for woven fabric
'weave type';                     fiberType = 'plainWeave';
'lamina thickness [mm]';     h_p = ht/nlayers; % 
'thickness of the fill [mm]';      h_f = h_p/2;
'thickness of the warp [mm]';      h_w = h_p/2;
'width of the fill [mm]';          a_f = 1.92;
'width of the warp [mm]';      a_w = 2;
'width of the fill gap [mm]';      g_f = 0.05;
'width of the warp gap [mm]';  g_w = 0.05;   
%%
% known parameters
m=8.55; % total mass of the specimen [kg]
V=1.2*1.2*ht; % specimen volume [m^3]
rho = m/V;
%% input for optimization
% lower and upper bounds of variables
%run(['inputs',filesep,'input1.m']);  % initial material properties
rhof0 = 1900; % kg/m^3
rhom0 = 1250; % kg/m^3
em0 = 3.43e9; % Pa
ef0 = 240e9; % Pa
nim0 = 0.35;
nif0 =  0.2; 
vol0 = 0.5;

e11_m = em0/1e9;
e11_f = ef0/1e9;
ni12_m = nim0;
ni12_f = nif0;
vol_0 = vol0;
rho_f = rhof0;
rho_m = rhom0;
e22_f = 0.1*e11_f;
ni23_f = ni12_f;
%% Mechanical properties  
 [C11,C12,C13,C21,C22,C23,C31,C32,C33,C44,C55,C66,rho] = ...
        compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_0, ...
        e11_m, ni12_m, rho_m, e11_f, e22_f, ni12_f, ni23_f, rho_f,false);
%% compute dispersion curves
[wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
%% compute objective function score
[score] = objective_fun(Data_polar,fmax,FREQ,number_of_modes_considered);
% fittnes function scaling factors
a=100;
b=310;   
ObjVal=a*(-1)*score+b;
if(~radians)
   wavenumber = wavenumber/(2*pi); % linear scale [1/m]
   wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
end
for j=1:number_of_angles % beta
    
    figfilename = [modelname,'_','angle_',num2str(beta(j)),'_dispersion_curves_initial_large'];
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
        fprintf('Making figure: [%d/%d]\n', j,number_of_angles);
        figure(j)
        set(gcf,'Color','w');
        %imagesc(fvec(2:end)/1e3, wavenumber(2:end,j)/2, squeeze(abs(Data_polar(j,2:end,2:end)))); 
        imagesc(fvec(2:end)/1e3, wavenumber(2:end,j), squeeze(abs(Data_polar(j,2:end,2:end)))); 
        set(gca,'YDir','normal'); 
        %axis([0 600 0 2000]);
        %axis([0 350 0 2000]);
        axis([0 350 0 min(wavenumber_max)]);
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
        fvec1=squeeze(FREQ(1,:,j)); % mode 1, angle j
        fvec2=squeeze(FREQ(2,:,j)); % mode 2, angle j
        fvec3=squeeze(FREQ(3,:,j)); % mode 3, angle j
        fvec4=squeeze(FREQ(4,:,j)); % mode 4, angle j
        kvec=squeeze(wavenumber(:,j)); % angle j
        %LW=0.5; % small figures
        LW=1; % large figures
        plot(fvec1(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
        plot(fvec2(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
        plot(fvec3(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
        plot(fvec4(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
        % get(gca,'FontName'); % default 'Helvetica'
        %set(gca,'FontName','Arial');
        %set(gca,'FontName','Helvetica');
        set(gca,'FontName','Times');
        fig = gcf;
        title({['$F=$ ',num2str(ObjVal),', ',num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
        print([output_path,figfilename],'-dpng', '-r600'); 
        %% END PLOTTING
    else
        fprintf('Figure: %s already exist\n', figfilename);
    end
close all;
end
%% END PLOTTING

