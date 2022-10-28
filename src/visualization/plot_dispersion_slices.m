%% Plot dispersion slices at 0, 90 deg and selected frequency


clear all;close all;   warning off;clc;
tic
load project_paths projectroot src_path;
set(0,'defaultAxesFontName', 'Times');
%%
% paths to full field measurements on NAS
list = {'Laminat_7_GFRP_Tomka/chirp_interp','aidd/data/raw/exp/L3_S4_B/chirp_interp','lamb_opt/data/raw/exp/chirp_interp'};
%list = {'Laminat_7_GFRP_Tomka/chirp_interp','lamb_opt/data/raw/exp/plain_weave/chirp_interp','lamb_opt/data/raw/exp/chirp_interp'};

% SANI - slightly anisotropic
% HANI - highly anisotropic
list_names={'SANI_GFRP','SANI_CFRP','HANI_CFRP'};
freq_list =[300,300,300]; % frequency list in kHz according to files above (max 4 frequencies)
%freq_list =[50,50,50]; % frequency list in kHz according to files above (max 4 frequencies)
test_case=[2,3]; % select file numbers for processing (starting from 2, chirp should be excluded)
%% Prepare output directories
% allow overwriting existing results if true
%overwrite=false;
overwrite=true;

currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
% prepare figure paths
figure_output_path = prepare_figure_paths(modelname);

radians_flag = false; % if true units of wanumbers [rad/m] if false [1/m]

%% input for figures
Cmap = jet(256); 
Cmap2 = turbo; 
caxis_cut = 0.8;
fig_width =6; % figure widht in cm
fig_height=5; % figure height in cm


%% Processing parameters
Nx = 512;   % number of points after interpolation in X direction
Ny = 512;   % number of points after interpolation in Y direction
N = 1024;% for zero padding

%%
% create path to the experimental raw data folder

raw_data_path = ['/pkudela_odroid_laser/'];

for k=test_case
    % load chirp signal
    disp('loading chirp signal');
    filename = list{k};
    load([raw_data_path,filename]); % Data, time, WL
    %Data=rot90(Data);
    [m,n,nft]=size(Data);
    Width=WL(1);
    Length=WL(2);

    disp('Transform to wavenumber-wavenumber-frequency domain');
    [KXKYF,kx_vec,ky_vec,f_vec] = spatial_to_wavenumber_wavefield_full2(Data,Length,Width,time); % full size data (-kx:+kx,-ky:+ky,-f:+f)
    [m1,n1,nft1] = size(KXKYF);
    % filter 3D wavefield for mode separation (A0 mode extraction)
    [kx_grid,ky_grid]=ndgrid(kx_vec,ky_vec);
    %figure;
    %surf(kx_grid,ky_grid,squeeze(abs(KXKYF(m+1:end,n+1:end,100))));shading interp; view(2);axis square
    [f_grid,k_grid]=ndgrid(f_vec,kx_vec);
    %figure; surf(f_grid,k_grid,squeeze(abs(KXKYF(m+1:end,m+1,nft+1:end)))');shading interp; view(2)
    %figure; surf(squeeze(abs(KXKYF(:,m+1,nft+1:end))));shading interp; view(2)
    fmax=f_vec(end);
    kxmax=kx_vec(end);
    kymax=ky_vec(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% KX-KY-F slices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(radians_flag)   
        [mkx,mky,mf] = meshgrid(kx_vec,ky_vec,f_vec/1e3);
    else
        [mkx,mky,mf] = meshgrid(kx_vec/(2*pi),ky_vec/(2*pi),f_vec/1e3);
    end
    % maxkx = 1000/(2*pi);
    % maxky = 1000/(2*pi);
    maxkx = 400;
    maxky = 400;
    maxf = 500;
    freq_slice = freq_list(k); % [kHz]
    xslice1 = []; yslice1 = []; zslice1 = freq_slice;
    xslice2 = 0; yslice2 = 0; zslice2 = [];
    figure;
    t=tiledlayout(2,1);
    %t.TileSpacing = 'tight';
    t.TileSpacing = 'none';
    t.Padding = 'tight';
    % Top plot
    ax1 = nexttile;
    h1 = slice(ax1,mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xslice2,yslice2,zslice2);
    set(h1,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');
    hold on;
    
    ylabel({'$k_y$ [1/m]'},'Rotation',-38,'Fontsize',12,'interpreter','latex');% for 7.5cm figure
    xlabel({'$k_x$ [1/m]'},'Rotation', 15,'Fontsize',12,'interpreter','latex');
    zlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex')
    set(ax1,'Fontsize',10,'FontName','Times','linewidth',1);
  
    grid(ax1,'off');
    view(3);
    lightangle(ax1,-45,45)
    lightangle(ax1,-45,45)
    colormap(Cmap2);
    %colormap turbo;
    line([0,0],[0,0],[0,max(f_vec)],'Color','y','LineWidth',1);
    line([-maxkx -maxkx],[-maxky  maxky],[freq_slice,freq_slice],'Color','r','LineWidth',1,'LineStyle','--');
    line([-maxkx  maxkx],[ maxky  maxky],[freq_slice,freq_slice],'Color','r','LineWidth',1,'LineStyle','--');
    line([ maxkx  maxkx],[ maxky -maxky],[freq_slice,freq_slice],'Color','r','LineWidth',1,'LineStyle','--');
    line([ maxkx -maxkx],[-maxky -maxky],[freq_slice,freq_slice],'Color','r','LineWidth',1,'LineStyle','--');
    xlim([-maxkx maxkx])
    ylim([-maxky maxky])
    zlim([0 maxf])
    box on; ax1.BoxStyle = 'full';
    %view(-20,20)
    %view(-40,15)
    view(-30,50)
    Smax=max(max(max(abs(KXKYF(3:end,3:end,end/2+10:end)))));
    %caxis([0 0.5*Smax]);
    %caxis([0 0.4*Smax]);
    caxis([0 0.3*Smax]);

    % bottom plot
    ax2 = nexttile;
    h2 = slice(ax2,mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xslice1,yslice1,zslice1);
    set(h2,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');
    hold on;
   
    ylabel({'$k_y$ [1/m]'},'Rotation',-38,'Fontsize',12,'interpreter','latex');% for 7.5cm figure
    xlabel({'$k_x$ [1/m]'},'Rotation', 15,'Fontsize',12,'interpreter','latex');
    zlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex')
    set(ax2,'Fontsize',10,'FontName','Times','linewidth',1);
    
    view(3);

    colormap(Cmap2);
    %colormap turbo;
%     lightangle(ax2,-45,45)
%     lightangle(ax2,-45,45)
    line([0,0],[0,0],[0,max(f_vec)],'Color','y','LineWidth',1);
    line([-maxkx -maxkx],[-maxky  maxky],[freq_slice,freq_slice],'Color','r','LineWidth',1,'LineStyle','--');
    line([-maxkx  maxkx],[ maxky  maxky],[freq_slice,freq_slice],'Color','r','LineWidth',1,'LineStyle','--');
    line([ maxkx  maxkx],[ maxky -maxky],[freq_slice,freq_slice],'Color','r','LineWidth',1,'LineStyle','--');
    line([ maxkx -maxkx],[-maxky -maxky],[freq_slice,freq_slice],'Color','r','LineWidth',1,'LineStyle','--');
    xlim([-maxkx maxkx])
    ylim([-maxky maxky])
    zlim([freq_slice-0.01*freq_slice freq_slice+0.01*freq_slice]);
    % box off; ax = gca; ax.BoxStyle = 'full';
    %axis on;
    axis off;
    grid(ax2,'off');
    %title([num2str(freq_slice),' kHz'],'Fontsize',10,'interpreter','latex');
    text(-maxkx,maxky,freq_slice+0.01*freq_slice,[num2str(freq_slice),' kHz'],'HorizontalAlignment','left','Fontsize',12,'interpreter','latex');
    %view(-20,20)
    %view(-40,50)
    view(-30,50)
    [~,I]=min(abs(freq_slice-f_vec/1e3));
    Smax=max(max(max(abs(KXKYF(3:end,3:end,end/2+I)))));
    %caxis([0 0.8*Smax]);
    %caxis([0 0.7*Smax]);
    %caxis([0 0.6*Smax]);
    caxis([0 0.3*Smax]);
    set(gcf,'color','white');set(gca,'TickDir','out');
    %set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
    set(gca, 'OuterPosition',[0 0 1. 1.]); % figure without axis and white border
    %set(gcf, 'Units','centimeters', 'Position',[10 10 8 10]);
    set(gcf, 'Units','centimeters', 'Position',[10 10 7.5 10]);
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));

    set(gcf,'PaperPositionMode','auto');
    drawnow;
    processed_filename = [list_names{k},'_',num2str(freq_slice),'_kHz']; % filename of processed .mat data
    print([figure_output_path,processed_filename],'-dpng', '-r600'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  