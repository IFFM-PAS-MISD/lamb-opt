% plot numerical dispersion curves

clear all; close all;

%set(0,'defaulttextinterpreter','none');

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true

% figure parameters
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
%fig_width = 7; fig_height = 5; 
% create path to the numerical model data folder
modelfolder = 'SASE';
modelname = 'SASE2';
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

for j=3 % beta = 30 deg
    
    figfilename = [modelname,'_','angle_',num2str(beta(j)),'_exemplary_dispersion_curves_sorted'];
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
        
        figure(j)
        set(gcf,'Color','w');
        for test_case=[7]
            % load numerical data file
            input_name = [model_input_path,num2str(test_case),'output'];
            load(input_name); % FREQ CG wavenumber
            wavenumber_step = wavenumber(3,j)-wavenumber(2,j);
            % automatic mode tracing - does not work
            %[cg_new,om_new] = mode_tracing_pade(CG(:,:,j),squeeze(FREQ(:,:,j)),wavenumber_step);
            % manual mode tracing
            [m,n,b] = size(FREQ);
            om_new=zeros(m,n);
            om_new(1,:) = FREQ(1,:,j);
            om_new(2,1:200) = FREQ(2,1:200,j);om_new(2,201:n) = FREQ(2,201:n,j);
            om_new(3,201:n) = FREQ(3,201:n,j);om_new(3,1:88) = FREQ(3,1:88,j);om_new(3,89:168) = FREQ(4,89:168,j); om_new(3,169:200) = FREQ(3,169:200,j); 
            om_new(4,1:88) = FREQ(4,1:88,j);om_new(4,89:168) = FREQ(3,89:168,j); om_new(4,169:n) = FREQ(4,169:n,j);
            if(~radians)
               wavenumber = wavenumber/(2*pi); % linear scale [1/m]
            end
          
            hold on;
%             fvec1=squeeze(FREQ(1,:,j)); % mode 1, angle j
%             fvec2=squeeze(FREQ(2,:,j)); % mode 2, angle j
%             fvec3=squeeze(FREQ(3,:,j)); % mode 3, angle j
%             fvec4=squeeze(FREQ(4,:,j)); % mode 4, angle j
            fvec1=squeeze(om_new(1,:)); % mode 1, angle j
            fvec2=squeeze(om_new(2,:)); % mode 2, angle j
            fvec3=squeeze(om_new(3,:)); % mode 3, angle j
            fvec4=squeeze(om_new(4,:)); % mode 4, angle j
            kvec=squeeze(wavenumber(:,j)); % angle j
            plot(fvec1(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[0,0,0]);
            plot(fvec2(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[1,0,0]);
            plot(fvec3(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[0,1,0]);
            plot(fvec4(2:end)/1e3,kvec(2:end),'linewidth',1,'color',[0,0,1]);
            
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
        title({[num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
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
%close all;
end
%% END PLOTTING

