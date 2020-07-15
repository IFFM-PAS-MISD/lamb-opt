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
fig_width = 7; fig_height = 7; 
% create path to the numerical model data folder
modelfolder = 'SASE';
modelname = 'SASE36_uni_surf';
radians = false;
% create output path
output_path = prepare_figure_paths(modelfolder,modelname);

% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname,'_out'], filesep );
beta=linspace(0,90,512);
selected_frequencies=[100,200,300,400]*1e3; % [Hz]
number_of_frequency_points=length(selected_frequencies);
% if(~radians)
%    wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
% end
fprintf('Making figures: %s\n', modelname);
% adjust colors
alpha_r = [0,0,0,0,0,0,1,1,1,1,1];
alpha_g = [0,0,0,0,0,0,0,0,0,0,0];
alpha_b = [1,1,1,1,1,0,0,0,0,0,0];
for j=1:number_of_frequency_points  % selected frequencies
    
    figfilename = [modelname,'_','frequency_',num2str(j),'_param_dispersion_surf_color'];
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
         fprintf('Making figure: [%d/%d]\n', j,number_of_frequency_points);
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
            dotsize = 3;
            plot(wx1(2:end),wy1(2:end),'.','color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)],'MarkerSize',dotsize);
            plot(wx2(2:end),wy2(2:end),'.','color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)],'MarkerSize',dotsize);
            plot(wx3(2:end),wy3(2:end),'.','color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)],'MarkerSize',dotsize);
            plot(wx4(2:end),wy4(2:end),'.','color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)],'MarkerSize',dotsize);
            plot(wx5(2:end),wy5(2:end),'.','color',[alpha_r(test_case),alpha_g(test_case),alpha_b(test_case)],'MarkerSize',dotsize);
        end
        box on;
        
        set(gca,'Fontsize',10,'linewidth',1);
        
        axis equal;
        axis([0 300 0 300]);
        set(gca,'Fontsize',10,'linewidth',1);
        if(radians)
            xlabel({'$k_x$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
            ylabel({'$k_y$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
        else
            xlabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
            ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
        end
        set(gca,'FontName','Times');
        fig = gcf;
        title({['${C_{44}} \pm$','30{\%}, ',num2str(selected_frequencies(j)/1000,'%5.0f'),' [kHz]']},'Fontsize',12,'interpreter','latex');
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

