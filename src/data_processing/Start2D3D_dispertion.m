clear all; clc; close all
%c = 1;
for c = 1:19
    switch c
        case 1
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\493x493p_0-350kHz_CHP160_x3_18Vpp_50Hz.mat')
            file = '493x493p_0-350kHz_CHP160_x3_18Vpp_50Hz.mat';    
        case 2
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_CHP160_x3_18Vpp_50Hz.mat')
            file = '289x289p_0-350kHz_CHP160_x3_18Vpp_50Hz.mat';
        case 3
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_CHP160_x3_18Vpp_500Hz.mat')
            file = '289x289p_0-350kHz_CHP160_x3_18Vpp_500Hz.mat';
        case 4
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_CHP80_x3_18Vpp_50Hz.mat')
            file = '289x289p_0-350kHz_CHP80_x3_18Vpp_50Hz.mat';
        case 5
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_CHP40_x3_18Vpp_50Hz.mat')
            file = '289x289p_0-350kHz_CHP40_x3_18Vpp_50Hz.mat';
        case 6
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_CHP20_x3_18Vpp_50Hz.mat')
            file = '289x289p_0-350kHz_CHP20_x3_18Vpp_50Hz.mat';
        case 7
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_CHP160_x3_18Vpp_5000Hz.mat')
            file = '289x289p_0-350kHz_CHP160_x3_18Vpp_5000Hz.mat';    %problem z PZT za male amplitudy        
        case 8
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_CHP160_x3_2Vpp_200Hz.mat')
            file = '289x289p_0-350kHz_CHP160_x3_2Vpp_200Hz.mat';    
        case 9
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_CHP160_x3_18Vpp_200Hz.mat')
            file = '289x289p_0-350kHz_CHP160_x3_18Vpp_200Hz.mat';   
        case 10
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_0-350kHz_HANN500_x3_18Vpp_50Hz.mat')
            file = '289x289p_HANN500_x3_18Vpp_50Hz.mat';   
        case 11
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_HANN100_x3_10Vpp_200Hz.mat')
            file = '289x289p_HANN100_x3_10Vpp_200Hz.mat';    
        case 12
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_10HC400kHz_x3_7Vpp_50Hz.mat')
            file = '289x289p_10HC400kHz_x3_7Vpp_50Hz.mat'; 
        case 13
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_HANN25_x3_10Vpp_200Hz.mat')
            file = '289x289p_HANN25_x3_10Vpp_200Hz.mat';    
        case 14
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_HANN50_x3_10Vpp_200Hz.mat')
            file = '289x289p_HANN50_x3_10Vpp_200Hz.mat'; 
        case 15
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_HANN50_x3_5Vpp_200Hz.mat')
            file = '289x289p_HANN50_x3_5Vpp_200Hz.mat';              
         case 16
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_HANN50_x3_20Vpp_200Hz.mat')
            file = '289x289p_HANN50_x3_15Vpp_200Hz.mat';            
         case 17
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_HANN50_x3_5Vpp_200Hz.mat')
            file = '289x289p_HANN50_x3_20Vpp_200Hz.mat';                 
         case 18
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\289x289p_HANN100_x30_10Vpp_200Hz.mat')
            file = '289x289p_HANN100_x30_10Vpp_200Hz.mat';    
         case 19
            point = 202; 
            load('l:\Praca\Pomiary Propagacja\Dispertion Curves CFRP\493x493p_HANN100_x10_10Vpp_200Hz.mat')
            file = '493x493p_HANN100_x10_10Vpp_200Hz.mat';              
    end
    


%Data(:,:,129:end) = [];
%time(1:2:end) = [];
%file(1:end-13) = [];

Lenght = 0.493;    
Width = 0.493;
cmap = jet(255);
[nY nX nT] = size(Data);

nsY = 1024;
nsX = 1024;
nsT = 1024;
siz = [nsY nsX nsT];


%% axis vectors
Fs =  1/(time(3)-time(2));                % sampling frequency
f_vec = Fs/2*linspace(0,1,nsT/2);         % frequency vector

dx = Lenght/(nX-1);
dkx = 1/(nX*dx);
kxmax = 1/(2*dx)-dkx/2;
kx_vec = 2*pi*linspace(0,kxmax,nsX/2);    % rad/m

dy = Width/(nY-1);
dky = 1/(nY*dy);
kymax = 1/(2*dy)-dky/2;
ky_vec = 2*pi*linspace(0,kymax,nsY/2);    % rad/m


%{
deltay = Lenght/(nY-1);       
deltax = Width/(nX-1);

ky1 = (0:nY-1)/nY*2*pi/deltay;              
ky_vec = (ky1-1/2*ky1(end));              % wavenumber rad/m
kx1 = (0:nX-1)/nX*2*pi/deltax;              
kx_vec = (kx1-1/2*kx1(end));              % wavenumber rad/m

y_vec = linspace(0,Lenght,nY);
x_vec = linspace(0,Width,nX);

clear ky1 kx1
%}



%% 3D FFT
KXKYF = fftshift(fftn(Data,siz));
%plot(abs(KXKYF(:,513,200)))

%% KY-freq
figure
KYF = squeeze(abs(KXKYF(floor(end/2)+1:end,floor(end/2)+1,floor(end/2)+1:end)));

%set(gcf,'Color','w','Position',[100 100 640 480]); 
imagesc(f_vec/1000, ky_vec/2/pi, KYF); set(gca,'YDir','normal');  

KYFrow = reshape(KYF,(floor(nsY/2))*(floor(nsT/2)),1);
KYFmax = max(KYFrow);
x_hist = 0:KYFmax/10000:KYFmax;
[n_elements xout] = hist(KYFrow,x_hist);
c_elements = cumsum(n_elements)/((floor(nsY/2))*(floor(nsT/2)));

i = 1;
while c_elements(i) < 0.999
    i = i+1;
end

colormap(cmap); axis tight
set(gca,'Fontsize',16)
b = max(caxis);  a = min(caxis); %b = 300;  
%caxis([xout(10) xout(i)]); 
caxis([0 xout(i)]); 

ylabel('k_y [1/m]'); xlabel('Frequency [kHz]'); 
title('3D FFT @ K_x = 0') 
%xlim([0 f_vec(end)/1000]); 
xlim([0 400]); 
%ylim([0 max(kx_vec/1000)])
ylim([0 250]);
set(gca,'FontName','Times');
fig = gcf;
width = 2*8;
height = 2*8;
set(fig, 'Units','centimeters', 'Position',[10 10 width height]); % size 12cm by 8cm (1-column text)
fig.PaperPositionMode  = 'auto';
set(gca,'LooseInset', max(get(gca,'TightInset'), 0));
print('-dpng','-r600',[file(1:end-4),'_KY_freq.png']);


%% KX-freq
%{
figure
set(gcf,'Color','w','Position',[100 100 640 480]); 
KXF =  abs(squeeze(abs(KXKYF(floor(end/2)+1,:,floor(end/2)+1:end))));
imagesc(f_vec/1000, kx_vec/1000, KXF); set(gca,'YDir','normal');  
colormap(cmap); axis tight
set(gca,'Fontsize',16)
ylabel('k_x (1/mm)'); xlabel('Frequency [kHz]'); 
title('3D FFT @ ky = 0') 
b = max(caxis);  
b = 300;
a = min(caxis);  caxis([0 b]); 
%xlim([0 f_vec(end)/1000]); 
xlim([0 350]); 
%ylim([0 max(kx_vec/1000)])
ylim([0 3.2]);
print('-djpeg','-r300',[file(1:end-4),'_KX2_freq.jpg']);
%}

%% KX-KY
%{
    point = 150;
    figure
for point = 150:150
    set(gcf,'Color','w','Position',[100 100 640 480]);
    imagesc(kx_vec/1000, ky_vec/1000, abs(KXKYF(:,:,end/2+1+point))); set(gca,'YDir','normal');  
    colormap(cmap); axis tight
    set(gca,'Fontsize',16)
    xlabel('k_x (1/mm)'); ylabel('k_y (1/mm)'); 
    title(['3D FFT @ ',num2str(round(f_vec(point)/1000)),' kHz']) 
    b = max(caxis);  a = min(caxis);  caxis([0 b/1.2]);  daspect([1,1,1])
    print('-djpeg','-r300',[file(1:end-4),'_KX_KY_',num2str(point),'.jpg']);
    %draw
    %wait(0.05)
end
%}

KXKYF_ = KXKYF(floor(end/2)+1:end,floor(end/2)+1:end,floor(end/2)+1:end); 
save([file(1:end-4), '_KXKYF'],'KXKYF_','f_vec','kx_vec', 'ky_vec');

end

%plot(squeeze(abs(KXKYF(50,50,:))))