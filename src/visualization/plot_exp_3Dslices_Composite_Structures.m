clc;clear;close all
load project_paths projectroot src_path;
%% Prepare output directories
% allow overwriting existing results if true
overwrite=false;
%overwrite=true;
% retrieve model name based on running file and folder
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
% idx = strfind( pathstr,filesep );
% modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 

tic
filename = '499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni.mat';

WL = [0.455 0.455];             % lengh width inb [m]
maxf = 500;                     % limit of the frequency in [kHz]
maxkx = 400;                    % limit Kx in [1/m]
maxky = 400;                    % limit Ky in [1/m]
%
freq_slice_list = [100 350]; 
angle_slice_list = [0 -60];% in respect to y axis
% create path to the experimental raw data folder
raw_data_path = fullfile( projectroot, 'data','raw','exp', filesep );
load([raw_data_path,filename]);
%fig_width = 12; fig_height = 8; 
fig_width = 9; fig_height = 6; 
modelfolder = 'Composite_Structures_uni';
% create output path
output_path = prepare_figure_paths(modelfolder,modelname);

%% x,y,t,kx,ky,f vectors and grids for plots
[rows, cols, samples] = size(Data);
x = linspace(0,WL(1),cols);     % [m]
y = linspace(0,WL(2),rows);     % [m]
[mx,my,mt] = meshgrid(x,y,time);

deltaR =WL(1)/(rows-1);
deltaC = WL(2)/(cols-1);

KX1 = (0:cols-1)/(cols)/deltaC;  
KX = (KX1-1/2*KX1(end));        % wavenumber centering [1/m]

KY1 = (0:rows-1)/(rows)/deltaR;  
KY = (KY1-1/2*KY1(end));        % wavenumber centering [1/m]

T = time(2)-time(1);
f = linspace(0,1/T/2,samples/2)*10^-3;

[mkx,mky,mf] = meshgrid(KX,KY,f);

%% X-Y-t slices
fig1 = figure; 
xslice = WL(1)/2; yslice = WL(2)/2;  zslice = [1*10^-4 7*10^-4];    %slices positions
h = slice(mx,my,mt,Data,xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(fig1,'Renderer','zbuffer');

xlabel({'$x$ [m]'},'Fontsize',12,'interpreter','latex');
ylabel({'$y$ [m]'},'Fontsize',12,'interpreter','latex');
zlabel({'$t$ [s]'},'Fontsize',12,'interpreter','latex');

lightangle(-45,45)
lightangle(-45,45)
colormap (jet(255))
caxis([-5 5]*10^-3)
box on; ax = gca; ax.BoxStyle = 'full';
view(-38.5,16)
axis tight

set(fig1, 'Units','centimeters', 'Position',[2 2 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig1.PaperPositionMode   = 'auto';
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
figfilename = ['X_Y_t_slice_',filename(1:end-4)];
%print([output_path,figfilename],'-dpng', '-r600'); 


fig2 = figure; 
xslice = 0.25; yslice = 0.25; zslice = [1*10^-4 7*10^-4];
h = slice(mx,my,mt,abs(Data),xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');

xlabel({'$x$ [m]'},'Fontsize',12,'interpreter','latex');
ylabel({'$y$ [m]'},'Fontsize',12,'interpreter','latex');
zlabel({'$t$ [s]'},'Fontsize',12,'interpreter','latex');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');

lightangle(-45,45)
lightangle(-45,45)
colormap (jet(255))
caxis([0 5]*10^-3)
box on; ax = gca; ax.BoxStyle = 'full';
view(-38.5,16)
axis tight

set(fig2, 'Units','centimeters', 'Position',[2 2 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig2.PaperPositionMode   = 'auto';
figfilename = ['X_Y_t_slice2_',filename(1:end-4)];
%print([output_path,figfilename],'-dpng', '-r600'); 

%% FFT3D
KXKYF = fftshift(fftn(Data));

%%  KX-KY-F slices
fig3 = figure;

xslice = []; yslice = []; zslice = freq_slice_list;
h = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');

ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex');
xlabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex');
zlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex')
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');

lightangle(-45,45)
lightangle(-45,45)
colormap (jet(255))
caxis([0 20]*10)
xlim([-maxkx maxkx])
ylim([-maxky maxky])
zlim([0 maxf])
box on; ax = gca; ax.BoxStyle = 'full';
view(-20,20)


set(fig3, 'Units','centimeters', 'Position',[2 2 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig3.PaperPositionMode   = 'auto';
figfilename = ['kx-ky-f_slices__',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 


fig4 = figure;
xslice = []; yslice = []; zslice = freq_slice_list;
h = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');

ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex');
xlabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex');
zlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex')
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');

lightangle(-45,45)
lightangle(-45,45)
lightangle(-10,45)
colormap (jet(255))
caxis([0 20]*10)
box on; ax = gca; ax.BoxStyle = 'full';
view(70,20)
xlim([-maxkx maxkx])
ylim([-maxky maxky])
zlim([0 maxf])

set(fig4, 'Units','centimeters', 'Position',[2 2 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig3.PaperPositionMode   = 'auto';
figfilename = ['kx-ky-f_slices2__',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 

%% slicing At Arbitrary Angles
%You can also create slices that are oriented in arbitrary planes. To do this, 
%Create a slice surface in the domain of the volume (surf, linspace).
%Orient this surface with respect to the axes (rotate).
%Get the XData, YData, and ZData of the surface (get).
%Use this data to draw the slice plane within the volume.

fig5 = figure; 

npoints = 500;                             % number of points (in both direction) for slicing plane - almost not increasing calculation time
meshh = meshgrid(0:maxf/(npoints-1):maxf);   

for k=1:length(angle_slice_list)
    hsp = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
    rotate(hsp,[0,0,1],angle_slice_list(k),[0,0,0]);
    xd{k} = get(hsp,'XData');
    yd{k} = get(hsp,'YData');
    zd{k} = get(hsp,'ZData');
    
end
for k=1:length(angle_slice_list)
    hs = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd{k},yd{k},zd{k});
    set(hs,'FaceColor','interp','EdgeColor','none'); 
    hold on;
end    
% hsp2 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
% rotate(hsp2,[0,0,1],15,[0,0,0])
%  xd2 = get(hsp2,'XData');
%  yd2 = get(hsp2,'YData');
%  zd2 = get(hsp2,'ZData'); 
% hsp3 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
% rotate(hsp3,[0,0,1],30,[0,0,0])
%  xd3 = get(hsp3,'XData');
%  yd3 = get(hsp3,'YData');
%  zd3 = get(hsp3,'ZData'); 
% hsp4 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
% rotate(hsp4,[0,0,1],45,[0,0,0])
%  xd4 = get(hsp4,'XData');
%  yd4 = get(hsp4,'YData');
%  zd4 = get(hsp4,'ZData'); 
% hsp5 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
% rotate(hsp5,[0,0,1],60,[0,0,0])
%  xd5 = get(hsp5,'XData');
%  yd5 = get(hsp5,'YData');
%  zd5 = get(hsp5,'ZData'); 
%  
% hsp9 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
% rotate(hsp9,[0,0,1],90,[0,0,0])
%  xd9 = get(hsp9,'XData');
%  yd9 = get(hsp9,'YData');
%  zd9 = get(hsp9,'ZData');  
 
% h1 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd1,yd1,zd1);
% set(h1,'FaceColor','interp','EdgeColor','none'); 
% hold on
%h2 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd2,yd2,zd2);
%set(h2,'FaceColor','interp','EdgeColor','none'); 
% h3 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd3,yd3,zd3);
% set(h3,'FaceColor','interp','EdgeColor','none'); 
% h4 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd4,yd4,zd4);
% set(h4,'FaceColor','interp','EdgeColor','none'); 
% h5 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd5,yd5,zd5);
% set(h5,'FaceColor','interp','EdgeColor','none'); 
%h9 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd9,yd9,zd9);
%set(h9,'FaceColor','interp','EdgeColor','none');

%h6 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),[],[],50);
%set(h6,'FaceColor','interp','EdgeColor','none');
set(gcf,'Renderer','zbuffer');

xlabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex');
ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex');
zlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');

%title({['what','ever']},'Fontsize',12,'interpreter','latex');
colormap (jet(255));
caxis([0 20]*10);
box on; ax = gca; ax.BoxStyle = 'full';
view(-20,20);
%view(-4,32);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
xlim([-maxkx maxkx])
ylim([-maxky maxky])
zlim([0 maxf])
lightangle(45,45);
lightangle(-20,45)
lightangle(10,80);
lightangle(-20,20);

set(fig5, 'Units','centimeters', 'Position',[2 2 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
fig5.PaperPositionMode   = 'auto';
figfilename = ['KXKYF_slices_2',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 

return;

%%%%%
fig6 = figure;
imagesc(f,KX,squeeze(abs(KXKYF(:,(end+1)/2+1,end/2+1:end))))
ylabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex');
xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
axis tight;
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
caxis([10 500])
xlim([0 maxf])
ylim([0 400])
colormap (jet(255));
set(gca,'YDir','normal')
set(fig6, 'Units','centimeters', 'Position',[2 2 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
fig6.PaperPositionMode   = 'auto';
figfilename = ['KX-F_',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 


fig7=figure;
imagesc(f,KY,squeeze(abs(KXKYF((end+1)/2+1,:,end/2+1:end))))
ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex');
xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
axis tight;
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
caxis([10 500])
xlim([0 maxf])
ylim([0 400])
colormap (jet(255));
set(gca,'YDir','normal')
set(fig7, 'Units','centimeters', 'Position',[2 2 fig_width fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
fig7.PaperPositionMode   = 'auto';
figfilename = ['KY-F_',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 

fig8=figure;
quarter = squeeze(  abs(  KXKYF( floor(end/2):end, floor(end/2):end, 200  )  )  ) ;  %f(200) ~= 250kHz
imagesc(KX(floor(end/2):end),KY(floor(end/2):end),quarter);
ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex');
xlabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex');
axis tight;
axis equal
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
caxis([10 100])
xlim([0 maxkx])
ylim([0 maxky])
colormap (jet(255));
set(gca,'YDir','normal')
set(fig8, 'Units','centimeters', 'Position',[2 2 fig_height fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
fig8.PaperPositionMode   = 'auto';
figfilename = ['KX-KY_',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 


fig9=figure;
imagesc(x,y,abs(Data(:,:,150)));
ylabel({'$y$ [m]'},'Fontsize',12,'interpreter','latex');
xlabel({'$x$ [m]'},'Fontsize',12,'interpreter','latex');
axis equal;
axis tight;
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
caxis([0 0.02])
colormap (jet(255));
set(gca,'YDir','normal')
set(fig9, 'Units','centimeters', 'Position',[2 2 fig_height fig_height]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
fig9.PaperPositionMode   = 'auto';
figfilename = ['x-y_',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 