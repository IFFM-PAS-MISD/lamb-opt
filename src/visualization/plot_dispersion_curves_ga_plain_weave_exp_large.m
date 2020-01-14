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


%%
[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);
wavenumber_min = zeros(number_of_angles,1); % minimal wavenumbers [1/m]
wavenumber_step=zeros(number_of_angles,1);
for j=1:number_of_angles
    wavenumber_step(j)=(wavenumber_max(j)-wavenumber_min(j))/(number_of_wavenumber_points-1); % wavenumber step [1/m]
end
wavenumber = zeros(number_of_wavenumber_points,number_of_angles);
for k=1:number_of_wavenumber_points
    wavenumber(k,:) = wavenumber_min+(k-1)*wavenumber_step;
end

if(~radians)
   wavenumber = wavenumber/(2*pi); % linear scale [1/m]
   wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
end
for j=1:number_of_angles % beta
    
    figfilename = [modelname,'_','angle_',num2str(beta(j)),'_dispersion_curves_exp_large'];
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

       
        set(gca,'FontName','Times');
        fig = gcf;
        title({[num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
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

