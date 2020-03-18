clc;clear;close all
tic
run d:\GIT\lamb-opt\config\config_matlab.m

filename = '499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni.mat';
WL = [0.455 0.455]; %lengh width

% create path to the experimental raw data folder
raw_data_path = fullfile( projectroot, 'data','raw','exp', filesep );
load([raw_data_path,filename]);

fig_width = 12; fig_height = 8; 

modelfolder = 'Experimental';
modelname = '';
% create output path
output_path = prepare_figure_paths(modelfolder,modelname);

%%
[rows cols samples] = size(Data);
x = linspace(0,WL(1),cols);
y = linspace(0,WL(2),rows);
[mx,my,mt] = meshgrid(x,y,time);


%% X-Y-t slices
xslice = WL(1)/2; yslice = WL(2)/2; zslice = [1*10^-4 7*10^-4];
figure
h = slice(mx,my,mt,Data,xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');
caxis([-5 5]*10^-3)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Time (s)')
lightangle(-45,45)
lightangle(-45,45)
colormap (jet(255))
box on
view(-38.5,16)
axis tight

xslice = [0.25]; yslice = [0.25]; zslice = [1*10^-4 7*10^-4];
figure
h = slice(mx,my,mt,abs(Data),xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');
caxis([0 5]*10^-3)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Time (s)')
lightangle(-45,45)
lightangle(-45,45)
colormap (jet(255))
box on
view(-38.5,16)
axis tight
%camzoom(1.1)
%camproj perspective

%% X-Y-t/2 slices

[mx2,my2,mt2] = meshgrid(x,y,time(1:512));

xslice = [0.25]; yslice = [0.25]; zslice = [1*10^-4];

figure
h = slice(mx2,my2,mt2,abs(Data(:,:,1:512)),xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');
caxis([0 20]*10^-3)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Time (s)')
lightangle(-45,45)
lightangle(-45,45)
colormap (jet(255))
box on
view(-38.5,35)
axis tight


%% kx-ky-f axis vectors
deltaR =WL(1)/(rows-1);
deltaC = WL(2)/(cols-1);

KX1 = (0:cols-1)/(cols)/deltaC;  
KX = (KX1-1/2*KX1(end)); % wavenumber centering [1/m]

KY1 = (0:rows-1)/(rows)/deltaR;  
KY = (KY1-1/2*KY1(end)); % wavenumber centering [1/m]

T = time(2)-time(1);
f = linspace(0,1/T/2,samples/2)*10^-3;

KXKYF = fftshift(fftn(Data));

%%  KX-KY-F slices
[mkx,mky,mf] = meshgrid(KX,KY,f);

xslice = [0]; yslice = []; zslice = [100 450];
figure
h = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');
caxis([0 20]*10)
xlabel('k_x (1/m)')
ylabel('k_y (1/m)')
zlabel('Freqency (kHz)')
lightangle(-45,45)
lightangle(-45,45)
colormap (jet(255))
box on
view(-20,20)
axis tight

xslice = []; yslice = 0; zslice = [100 450];
figure
tic
h = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xslice,yslice,zslice);
set(h,'FaceColor','interp','EdgeColor','none'); set(gcf,'Renderer','zbuffer');
caxis([0 20]*10)
xlabel('k_x (rad/mm)')
ylabel('k_y (rad/mm)')
zlabel('Freqency (kHz)')
lightangle(-45,45)
lightangle(-45,45)
colormap (jet(255))
box on
view(70,20)
axis tight

%%
%datasmall = abs(KXKYF_(1:10:512,1:10:512,1:10:512));

%%
%slicing At Arbitrary Angles
%You can also create slices that are oriented in arbitrary planes. To do this, 
%Create a slice surface in the domain of the volume (surf, linspace).
%Orient this surface with respect to the axes (rotate).
%Get the XData, YData, and ZData of the surface (get).
%Use this data to draw the slice plane within the volume.
%For example, these statements slice the volume in the first example with a rotated plane. Placing these commands within a for loop "passes" the plane through the volume along the z-axis.

figure 
npoints = 500;                             %number of points (in both direction) for slicing plane - almost not increasing calculation time
maxf = 500; % limit of the frequency in kHz
meshh = meshgrid(0:maxf/(npoints-1):maxf);   

hsp1 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
rotate(hsp1,[0,0,1],0,[0,0,0])
 xd1 = get(hsp1,'XData');
 yd1 = get(hsp1,'YData');
 zd1 = get(hsp1,'ZData');
hsp2 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
rotate(hsp2,[0,0,1],15,[0,0,0])
 xd2 = get(hsp2,'XData');
 yd2 = get(hsp2,'YData');
 zd2 = get(hsp2,'ZData'); 
hsp3 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
rotate(hsp3,[0,0,1],30,[0,0,0])
 xd3 = get(hsp3,'XData');
 yd3 = get(hsp3,'YData');
 zd3 = get(hsp3,'ZData'); 
hsp4 = surf(linspace(0,0,npoints),linspace(KY(1),KY(end),npoints),meshh);
rotate(hsp4,[0,0,1],45,[0,0,0])
 xd4 = get(hsp4,'XData');
 yd4 = get(hsp4,'YData');
 zd4 = get(hsp4,'ZData'); 
 
h1 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd1,yd1,zd1);
set(h1,'FaceColor','interp','EdgeColor','none'); 
hold on
%h2 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd2,yd2,zd2);
%set(h2,'FaceColor','interp','EdgeColor','none'); 
%h3 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd3,yd3,zd3);
%set(h3,'FaceColor','interp','EdgeColor','none'); 
h4 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),xd4,yd4,zd4);
set(h4,'FaceColor','interp','EdgeColor','none'); 
h5 = slice(mkx,mky,mf,abs(KXKYF(:,:,end/2+1:end)),[],[],[50]);
set(h5,'FaceColor','interp','EdgeColor','none');
set(gcf,'Renderer','zbuffer');
caxis([0 20]*10);
xlabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex');
ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex');
zlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
%title({['what','ever']},'Fontsize',12,'interpreter','latex');
lightangle(-45,45);
lightangle(-45,45);
colormap (jet(255));
box on;
view(-20,20);
%view(-4,32);
axis tight;
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
zlim([0 maxf])

fig = gcf;
set(fig, 'Units','centimeters', 'Position',[2 2 fig_width*3 fig_height*3]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.1));
fig.PaperPositionMode   = 'auto';
figfilename = ['KXKYF_angles_slice_2small',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 

%%%%%
figure
imagesc(f,KX,squeeze(abs(KXKYF(:,(end+1)/2,end/2+1:end))))
xlabel({'$k_x$ [1/m]'},'Fontsize',12,'interpreter','latex');
ylabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
axis tight;
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
caxis([0 500])
xlim([0 500])
ylim([0 300])
set(gca,'YDir','normal')
fig = gcf;
set(fig, 'Units','centimeters', 'Position',[2 2 fig_width*3 fig_height*3]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
fig.PaperPositionMode   = 'auto';
figfilename = ['KYF',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 


figure
imagesc(f,KY,squeeze(abs(KXKYF((end+1)/2,:,end/2+1:end))))
xlabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex');
ylabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
axis tight;
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
caxis([0 500])
xlim([0 500])
ylim([0 300])
set(gca,'YDir','normal')
fig = gcf;
set(fig, 'Units','centimeters', 'Position',[2 2 fig_width*3 fig_height*3]); % size 12cm by 8cm (1-column text)
% remove unnecessary white space
fig.PaperPositionMode   = 'auto';
figfilename = ['KYF',filename(1:end-4)];
print([output_path,figfilename],'-dpng', '-r600'); 



