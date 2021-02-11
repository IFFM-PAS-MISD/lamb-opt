% plot numerical dispersion curves (kx-ky surface plot)

clear all; close all;

%set(0,'defaulttextinterpreter','none');

% load projectroot path
load project_paths projectroot src_path;
overwrite = true; % allow overwriting existing results if true

% figure parameters
% size 12cm by 12cm (1-column text)
%fig_width = 12; fig_height = 12; 
% size 7cm by 7cm (2-column text)
fig_width = 7; fig_height = 7; 
% create path to the numerical model data folder
modelfolder = 'genetic_algorithm';
modelname = 'ga_unidirectional_C_tensor_known_mass_kx_ky_stringer';
radians = false;
test_case=1; % numerical data
% create output path
output_path = prepare_figure_paths(modelfolder,modelname);

% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname,'_out'], filesep );
% create path to the experimental interim data folder
exp_input_data_path=fullfile( projectroot, 'data','interim','exp', filesep );
% create path to the experimental processed data folder
exp_param_data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
input_file = 1; % experimental data
% load experimental data file
exp_filename = {'interim_483x483p_CHIRP_20-500kHz_125us_6Vpp_x3_stringer_intact_KXKYF'};     

load([exp_input_data_path,exp_filename{input_file}]); % Data wavenumber_max fmax beta number_of_wavenumber_points  
% filename of parameter data
filename = {'polar_f_interim_483x483p_CHIRP_20-500kHz_125us_6Vpp_x3_stringer_intact_KXKYF_param'}; 

load([exp_param_data_path,filename{input_file}]); % beta number_of_wavenumber_points selected_frequency_index selected_frequencies ...
fmin = selected_frequencies(1);
fmax = selected_frequencies(end);
number_of_frequency_points = length(selected_frequencies);
%% Input for SASE

np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
%layup = [45,0,-45,90,90,-45,0,45]-90; ht = 2.086/1000; % [m] laminate total thickness;% 1
%layup = [45,0,-45,90,-45,0,45,90,90,45,0,-45,90,-45,0,45]-90; ht = 2.086/1000; % [m] laminate total thickness;% 2
%layup = [45,0,-45,90,90,-45,0,45]; ht = 2.086/1000; % 3
layup = [45,0,-45,90,-45,0,45,90,90,45,0,-45,90,-45,0,45]; ht = 2.0/1000; % [m] laminate total thickness;% 4

nlayers = length(layup);

h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers; new plain weave
% Stacking direction
stack_dir = 1;
%%
% known parameters
m=6.46; % total mass of the specimen [kg]
V=1.2*1.2*ht; % specimen volume [m^3]
rho = 1571; % [kg/m^3]
%%

Data = KXKYF_(:,:,selected_frequency_index);
clear KXKYF_;
[number_of_wavenumber_points_x,number_of_wavenumber_points_y,~] = size(Data);



% if(~radians)
%    wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
% end
fprintf('Making figures: %s\n', modelname);
%% load optimized constants
output_name = [model_input_path,filesep,num2str(test_case),'output'];
load(output_name); % 'C11','C12','C13','C22','C23','C33','C44','C55','C66','rho','ObjVal'

%% compute dispersion curves
[wavenumber,CG,FREQ] = main_SASE2(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,fmin,fmax,number_of_frequency_points,beta,stack_dir,np,nele_layer);
% dispersion surf for reference data from Jochen article
[wavenumber_ref,CG_ref,FREQ_ref] = main_SASE2(rho,130e9,6.06e9,6.06e9,11.19e9,5.19e9,11.19e9,3e9,4.13e9,4.13e9,layup,h,fmin,fmax,number_of_frequency_points,beta,stack_dir,np,nele_layer);

% save file - it needs tweaking - add checking if C tensor was already computed 
% load('W.mat');
if(~radians)
   wavenumber = wavenumber/(2*pi); % linear scale [1/m]
   wavenumber_ref = wavenumber_ref/(2*pi); % linear scale [1/m]
   wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
end

kvecx = linspace(0,wavenumber_max,number_of_wavenumber_points_x);
kvecy = linspace(0,wavenumber_max,number_of_wavenumber_points_y);
for j=1:number_of_frequency_points  % selected frequencies
    
    figfilename = [modelname,'_','frequency_',num2str(j),'_dispersion_surf_test_case_',num2str(test_case),'_small4'];
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
        fprintf('Making figure: [%d/%d]\n', j,number_of_frequency_points);
        figure(j)
        set(gcf,'Color','w');
        imagesc(kvecx(2:end), kvecx(2:end), squeeze(abs(Data(2:end,2:end,j)))); 
        set(gca,'YDir','normal'); 
        axis([0 300 0 300]);
        set(gca,'Fontsize',10,'linewidth',1);
        if(radians)
            xlabel({'$k_x$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
            ylabel({'$k_y$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
        else
            xlabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
            ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
        end
        colormap jet; 
        %axis tight; 
        caxis([0 max(caxis)/3]); 
        % transform numerical dispersion surface from (k-beta) polar coordinate system to (kx-ky) Cartesian coordinate system
        wx1 = squeeze(wavenumber(1 , j, :)).*cos(beta'*pi/180); % mode 1, frequency j
        wy1 =  squeeze(wavenumber(1, j, :)).*sin(beta'*pi/180); %  mode 1, frequency j
        wx2 = squeeze(wavenumber(2 , j, :)).*cos(beta'*pi/180); % mode 2, frequency j
        wy2 =  squeeze(wavenumber(2, j, :)).*sin(beta'*pi/180); %  mode 2, frequency j
        wx3 = squeeze(wavenumber(3 , j, :)).*cos(beta'*pi/180); % mode 3, frequency j
        wy3 =  squeeze(wavenumber(3, j, :)).*sin(beta'*pi/180); %  mode 3, frequency j
        wx4 = squeeze(wavenumber(4 , j, :)).*cos(beta'*pi/180); % mode 4, frequency j
        wy4 =  squeeze(wavenumber(4, j, :)).*sin(beta'*pi/180); %  mode 4, frequency j
        wx5 = squeeze(wavenumber(5 , j, :)).*cos(beta'*pi/180); % mode 5, frequency j
        wy5 =  squeeze(wavenumber(5, j, :)).*sin(beta'*pi/180); %  mode 5, frequency j
        
        wx1ref = squeeze(wavenumber_ref(1 , j, :)).*cos(beta'*pi/180); % mode 1, frequency j
        wy1ref =  squeeze(wavenumber_ref(1, j, :)).*sin(beta'*pi/180); %  mode 1, frequency j
        wx2ref = squeeze(wavenumber_ref(2 , j, :)).*cos(beta'*pi/180); % mode 2, frequency j
        wy2ref =  squeeze(wavenumber_ref(2, j, :)).*sin(beta'*pi/180); %  mode 2, frequency j
        wx3ref = squeeze(wavenumber_ref(3 , j, :)).*cos(beta'*pi/180); % mode 3, frequency j
        wy3ref =  squeeze(wavenumber_ref(3, j, :)).*sin(beta'*pi/180); %  mode 3, frequency j
        wx4ref = squeeze(wavenumber_ref(4 , j, :)).*cos(beta'*pi/180); % mode 4, frequency j
        wy4ref =  squeeze(wavenumber_ref(4, j, :)).*sin(beta'*pi/180); %  mode 4, frequency j
        wx5ref = squeeze(wavenumber_ref(5 , j, :)).*cos(beta'*pi/180); % mode 5, frequency j
        wy5ref =  squeeze(wavenumber_ref(5, j, :)).*sin(beta'*pi/180); %  mode 5, frequency j
        hold on;
    
        LW=0.1; % small figures
        %LW=1; % large figures
        dotsize = 3;
        plot(wx1(2:end),wy1(2:end),'w.','MarkerSize',dotsize);
        plot(wx2(2:end),wy2(2:end),'w.','MarkerSize',dotsize);
        plot(wx3(2:end),wy3(2:end),'w.','MarkerSize',dotsize);
        plot(wx4(2:end),wy4(2:end),'w.','MarkerSize',dotsize);
        plot(wx5(2:end),wy5(2:end),'w.','MarkerSize',dotsize);
        % reference data from Jochen article
        plot(wx1ref(2:end),wy1ref(2:end),'m.','MarkerSize',dotsize);
        plot(wx2ref(2:end),wy2ref(2:end),'m.','MarkerSize',dotsize);
        plot(wx3ref(2:end),wy3ref(2:end),'m.','MarkerSize',dotsize);
        plot(wx4ref(2:end),wy4ref(2:end),'m.','MarkerSize',dotsize);
        plot(wx5ref(2:end),wy5ref(2:end),'m.','MarkerSize',dotsize);
        % get(gca,'FontName'); % default 'Helvetica'
        %set(gca,'FontName','Arial');
        %set(gca,'FontName','Helvetica');
        set(gca,'FontName','Times');
        fig = gcf;
        title({['$F=$ ',num2str(ObjVal,'%5.2f'),', ',num2str(selected_frequencies(j)/1000,'%5.2f'),' [kHz]']},'Fontsize',12,'interpreter','latex');
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

