% plot numerical dispersion curves

clear all; close all;

%set(0,'defaulttextinterpreter','none');

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true
%overwrite = true;
% figure parameters
% size 12cm by 8cm (1-column text)
%fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
fig_width = 7; fig_height = 5; 
% create path to the numerical model data folder
modelfolder = 'genetic_algorithm';
modelname = 'ga_unidirectional_C_tensor_known_mass_mut_rnd_offspring_2lay6';
modelname2= 'ga_uni_SHMII_2021_homogenized_num';
radians = false;
test_case=1; % numerical data
% create output path
output_path = prepare_figure_paths(modelfolder,modelname2);

% create path to the numerical raw data folder (ga)
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname,'_out'], filesep );
% create path to the numericall processed data folder (polar data)
num_input_data_path=fullfile( projectroot, 'data','processed','num', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
% and conversion to polar coordinate system

% filename of parameter data
filename = {'polar_interim_flat_shell_Vz_19_982x982bottom_KXKYF_param',...
                'polar_interim_flat_shell_Vz_20_982x982bottom_KXKYF_param',...
                'polar_interim_flat_shell_Vz_21_982x982bottom_KXKYF_param'}; 

% load numerical data file
num_filename = {'polar_interim_flat_shell_Vz_19_982x982bottom_KXKYF',...
                'polar_interim_flat_shell_Vz_20_982x982bottom_KXKYF',...
                'polar_interim_flat_shell_Vz_21_982x982bottom_KXKYF'}; 

input_file=1;
load([num_input_data_path,filename{input_file}]); % wavenumber_max fmax beta number_of_wavenumber_points           
load([num_input_data_path,num_filename{input_file}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points
wavenumber_max1 = wavenumber_max;
fmax1 =fmax;

Data_polar1=Data_polar;

input_file=2;
load([num_input_data_path,filename{input_file}]); % wavenumber_max fmax beta number_of_wavenumber_points           
load([num_input_data_path,num_filename{input_file}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points
wavenumber_max2 = wavenumber_max;
fmax2 =fmax;

Data_polar2=Data_polar;

input_file=3;
load([num_input_data_path,filename{input_file}]); % wavenumber_max fmax beta number_of_wavenumber_points           
load([num_input_data_path,num_filename{input_file}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points
wavenumber_max3 = wavenumber_max;
fmax3 =fmax;

Data_polar3=Data_polar;

clear Data_polar;

%% Input for SASE
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
ht = 2.85/1000; % [m] laminate total thickness; unidirectional
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
layup = [90 90 90 90 90 90 90 90];

nlayers = length(layup);

h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers; new plain weave
% Stacking direction
stack_dir = 1;
%%
[number_of_angles,number_of_wavenumber_points1,number_of_frequency_points1] = size(Data_polar1);
fvec1a = linspace(0,fmax1,number_of_frequency_points1);
[number_of_angles,number_of_wavenumber_points2,number_of_frequency_points2] = size(Data_polar2);
fvec2a = linspace(0,fmax2,number_of_frequency_points2);
[number_of_angles,number_of_wavenumber_points3,number_of_frequency_points3] = size(Data_polar3);
fvec3a = linspace(0,fmax3,number_of_frequency_points3);

fprintf('Making figures: %s\n', modelname);
%% load optimized constants
output_name = [model_input_path,filesep,num2str(test_case),'output'];
load(output_name); % 'rhom','rhof','em','ef','nim','nif','vol','C11','C12','C13','C22','C23','C33','C44','C55','C66','rho','ObjVal'
%% compute dispersion curves
[wavenumber1,CG1,FREQ1] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max1,number_of_wavenumber_points1,beta,stack_dir,np,nele_layer);
[wavenumber2,CG2,FREQ2] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max2,number_of_wavenumber_points2,beta,stack_dir,np,nele_layer);
[wavenumber3,CG3,FREQ3] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max3,number_of_wavenumber_points3,beta,stack_dir,np,nele_layer);

if(~radians)
   wavenumber1 = wavenumber1/(2*pi); % linear scale [1/m]
   wavenumber2 = wavenumber2/(2*pi); % linear scale [1/m]
   wavenumber3 = wavenumber3/(2*pi); % linear scale [1/m]
   wavenumber_max1 = wavenumber_max1/(2*pi); % linear scale [1/m]
   wavenumber_max2 = wavenumber_max2/(2*pi); % linear scale [1/m]
   wavenumber_max3 = wavenumber_max3/(2*pi); % linear scale [1/m]
end
[Xq,Yq]=meshgrid(linspace(0.1,200,800),linspace(0.1,120,480));

for j=1:number_of_angles % beta

    figfilename = [modelname2,'_','angle_',num2str(beta(j)),'_dispersion_curves_test_case_',num2str(test_case),'_small'];
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% compute combined dispersion image
        [X1,Y1] = meshgrid(wavenumber1(2:end,j),fvec1a(2:end)/1e3);
        V=squeeze( abs(Data_polar1(j,2:end,2:end)))';

        Vq1 = interp2(X1,Y1,V,Xq,Yq,'cubic',0);
        Vq1=Vq1/max(max(Vq1));


        [X2,Y2] = meshgrid(wavenumber2(2:end,j),fvec2a(2:end)/1e3);
        V=squeeze( abs(Data_polar2(j,2:end,2:end)))';
        Vq2 = interp2(X2,Y2,V,Xq,Yq,'cubic',0);
        Vq2=Vq2/max(max(Vq2));

        [X3,Y3] = meshgrid(wavenumber3(2:end,j),fvec3a(2:end)/1e3);
        V=squeeze( abs(Data_polar3(j,2:end,2:end)))';
        Vq3 = interp2(X3,Y3,V,Xq,Yq,'cubic',0);
        Vq3=Vq3/max(max(Vq3));
        Vq=(Vq1+Vq2+Vq3)/3; % sum of dispersion maps (wavepackets centered at 16.5, 50 and 100 kHz)
          
        %% START PLOTTING
        fprintf('Making figure: [%d/%d]\n', j,number_of_angles);
        figure(j)
        set(gcf,'Color','w');
     
        surf(Yq,Xq,Vq); shading interp; view(2);set(gcf, 'Renderer', 'zbuffer'); % numerical dispersion image
        
        shading interp; view(2);
        set(gcf, 'Renderer', 'zbuffer');
        axis([0 120 0 150]);
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
        fvec1=squeeze(FREQ3(1,:,j)); % mode 1, angle j
        fvec2=squeeze(FREQ3(2,:,j)); % mode 2, angle j
        fvec3=squeeze(FREQ3(3,:,j)); % mode 3, angle j
        fvec4=squeeze(FREQ3(4,:,j)); % mode 4, angle j
        fvec5=squeeze(FREQ3(5,:,j)); % mode 5, angle j
        kvec=squeeze(wavenumber3(:,j)); % angle j
        LW=0.5; % small figures
        %LW=1; % large figures
        plot3(fvec1(2:end)/1e3,kvec(2:end),zeros(length(kvec)-1,1)+1,'w','linewidth',LW); % A0 mode
%         plot3(fvec2(2:end)/1e3,kvec(2:end),zeros(length(kvec)-1,1)+1,'w','linewidth',LW);
%         plot3(fvec3(2:end)/1e3,kvec(2:end),zeros(length(kvec)-1,1)+1,'w','linewidth',LW);
%         plot3(fvec4(2:end)/1e3,kvec(2:end),zeros(length(kvec)-1,1)+1,'w','linewidth',LW);
%         plot3(fvec5(2:end)/1e3,kvec(2:end),zeros(length(kvec)-1,1)+1,'w','linewidth',LW);
        % get(gca,'FontName'); % default 'Helvetica'
        %set(gca,'FontName','Arial');
        %set(gca,'FontName','Helvetica');
        box on;
        set(gca,'FontName','Times');
        fig = gcf;
        title({['$F=$ ',num2str(ObjVal,'%5.2f'),', ',num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
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

