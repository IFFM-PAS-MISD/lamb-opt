% plot numerical dispersion curves on top of experimental data

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
foldername = 'minigrant_raport_2';
% create output path
output_path = prepare_exp_figure_paths(foldername);

% input and chirp signals
f0=0;
f1=350000;%instantaneous frequency f1 [Hz]
t1=20e-6; % at time t1 [s]
t20 = 0:t1/1e3:t1;
y20 = mychirp(t20,f0,t1,f1,pi/2);
t1=40e-6;
t40 = 0:t1/1e3:t1;
y40 = mychirp(t40,f0,t1,f1,pi/2);
t1=80e-6;
t80 = 0:t1/1e3:t1;
y80 = mychirp(t80,f0,t1,f1,pi/2);
t1=160e-6;
t160 = 0:t1/1e3:t1;
y160 = mychirp(t160,f0,t1,f1,pi/2);


fprintf('Making chirp figures for: %s\n', foldername); 
figfilename = 'chirp20';
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
        figure(1)  
        plot(t20*1e6,y20,'linewidth',1,'color','k');
        box on;
        axis([0 t20(end)*1e6 -1 1]);
        set(gca, 'Layer', 'Top');
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$t$ [$\mu$s]'},'Fontsize',12,'interpreter','latex');
        ylabel({'$A$ [-]'},'Fontsize',12,'interpreter','latex'); 
        set(gca,'FontName','Times');
        fig = gcf;
        
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]);
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
        print([output_path,figfilename],'-dpng', '-r600'); 
    else
        fprintf('Figure: %s already exist\n', figfilename);
    end
figfilename = 'chirp40';
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
        figure(2)  
        plot(t40*1e6,y40,'linewidth',1,'color','k');
        box on;
        axis([0 t40(end)*1e6 -1 1]);
        set(gca, 'Layer', 'Top');
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$t$ [$\mu$s]'},'Fontsize',12,'interpreter','latex');
        ylabel({'$A$ [-]'},'Fontsize',12,'interpreter','latex'); 
        set(gca,'FontName','Times');
        fig = gcf;
        
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]);
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
        print([output_path,figfilename],'-dpng', '-r600'); 
    else
        fprintf('Figure: %s already exist\n', figfilename);
    end
figfilename = 'chirp80';
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
        figure(2)  
        plot(t80*1e6,y80,'linewidth',1,'color','k');
        box on;
        axis([0 t80(end)*1e6 -1 1]);
        set(gca, 'Layer', 'Top');
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$t$ [$\mu$s]'},'Fontsize',12,'interpreter','latex');
        ylabel({'$A$ [-]'},'Fontsize',12,'interpreter','latex'); 
        set(gca,'FontName','Times');
        fig = gcf;
        
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]);
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
        print([output_path,figfilename],'-dpng', '-r600'); 
    else
        fprintf('Figure: %s already exist\n', figfilename);
    end
figfilename = 'chirp160';
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        %% START PLOTTING
        figure(2)  
        plot(t160*1e6,y160,'linewidth',1,'color','k');
        box on;
        axis([0 t160(end)*1e6 -1 1]);
        xticks([0 40 80 120 160]);
        set(gca, 'Layer', 'Top');
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$t$ [$\mu$s]'},'Fontsize',12,'interpreter','latex');
        ylabel({'$A$ [-]'},'Fontsize',12,'interpreter','latex'); 
        set(gca,'FontName','Times');
        fig = gcf;
        
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]);
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
        print([output_path,figfilename],'-dpng', '-r600'); 
    else
        fprintf('Figure: %s already exist\n', figfilename);
    end
close all;
%% END PLOTTING

