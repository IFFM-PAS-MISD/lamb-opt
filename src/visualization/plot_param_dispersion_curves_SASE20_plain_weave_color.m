% plot numerical dispersion curves

clear all; close all;

%set(0,'defaulttextinterpreter','none');

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true

% figure parameters
% size 12cm by 8cm (1-column text)
%fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
fig_width = 7; fig_height = 5; 
% create path to the numerical model data folder
modelfolder = 'SASE';
modelname = 'SASE20_plain_weave';
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
filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_param'; 

% load experimental parameter data file
load([exp_input_data_path,filename]); % wavenumber_max fmax beta number_of_wavenumber_points
number_of_frequency_points = number_of_wavenumber_points;
number_of_angles = length(beta);
fvec = linspace(0,fmax,number_of_frequency_points);
if(~radians)
   wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
end
fprintf('Making figures: %s\n', modelname);
% adjust colors
alpha_r = [0,0,0,0,0,0,1,1,1,1,1];
alpha_g = [0,0,0,0,0,0,0,0,0,0,0];
alpha_b = [1,1,1,1,1,0,0,0,0,0,0];
for j=1:number_of_angles % beta
    
    figfilename = [modelname,'_','angle_',num2str(beta(j)),'_param_dispersion_curves_color'];
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
        fprintf('Making figure: [%d/%d]\n', j,number_of_angles);
        figure(j)
        set(gcf,'Color','w');
        for test_case=[1:5,7:11,6]
            % load numerical data file
            input_name = [model_input_path,num2str(test_case),'output'];
            load(input_name); % FREQ CG wavenumber
            if(~radians)
               wavenumber = wavenumber/(2*pi); % linear scale [1/m]
            end
          
            hold on;
            fvec1=squeeze(FREQ(1,:,j)); % mode 1, angle j
            fvec2=squeeze(FREQ(2,:,j)); % mode 2, angle j
            fvec3=squeeze(FREQ(3,:,j)); % mode 3, angle j
            fvec4=squeeze(FREQ(4,:,j)); % mode 4, angle j
            kvec=squeeze(wavenumber(:,j)); % angle j
            plot(fvec1(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)]);
            plot(fvec2(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)]);
            plot(fvec3(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)]);
            plot(fvec4(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)]);
            
        end
        box on;
        %axis([0 350 0 min(wavenumber_max)]);
        axis([0 500 0 min(wavenumber_max)]);
        set(gca, 'Layer', 'Top');
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
        if(radians)
            ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
        else
            ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
        end
        set(gca,'FontName','Times');
        fig = gcf;
        title({['${C_{11}} \pm$','20{\%}, ',num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
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

