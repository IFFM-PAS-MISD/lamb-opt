% plot numerical dispersion curves (kx-ky surface plot)

clear all; close all;

%set(0,'defaulttextinterpreter','none');

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true

% figure parameters
% size 12cm by 12cm (1-column text)
%fig_width = 12; fig_height = 12; 
% size 7cm by 7cm (2-column text)
fig_width = 7; fig_height = 7; 
% create path to the numerical model data folder
modelfolder = 'genetic_algorithm';
modelname = 'ga_unidirectional_C_tensor_known_mass_kx_ky';
radians = false;
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder2 = pathstr(idx(end)+1:end); % name of folder
modelname2 = name; 

test_case=2; % numerical data
% path to dispersion curves data
dispersion_output_path = prepare_figure_paths(modelfolder,modelname);
% create output path
figure_output_path = prepare_figure_paths(modelname2);
output_path = prepare_folder_paths('processed','exp',modelname2);

% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname,'_out'], filesep );
% create path to the experimental interim data folder
exp_input_data_path=fullfile( projectroot, 'data','interim','exp', filesep );
% create path to the experimental processed data folder
exp_param_data_path=fullfile( projectroot, 'data','processed','exp', filesep );
% create path to the experimental wavefield data folder
exp_wavefield_data_path=fullfile( projectroot, 'data','raw','exp', filesep );
% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
input_file = 1; % experimental data
% disp('loading experimental data');
% % load experimental data file
% exp_filename = {'interim_499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni_KXKYF',... % 1 small area unidirectional
%                             'interim_499x499p_chp200_x40_6Vpp_250Hz_uni_KXKYF'};         % 2 large area unidirectional
% load([exp_input_data_path,exp_filename{input_file}]); % 'KXKYF_','f_vec','kx_vec', 'ky_vec' 

%% Input for SASE
ht = 2.85/1000; % [m] laminate total thickness; unidirectional
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
no_of_angles=512;
beta = 0:90/(no_of_angles-1):90;
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
wavenumber_max = zeros(length(beta),1)+3.43e3; % maximal wavenumbers [1/m]
number_of_wavenumber_points=512;
layup = [90 90 90 90 90 90 90 90];
nlayers = length(layup);
h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers; new plain weave
% Stacking direction
stack_dir = 1;
% known parameters
m=6.46; % total mass of the specimen [kg]
V=1.2*1.2*ht; % specimen volume [m^3]
rho = m/V;
%% input for mask
mask_width_A0_1=60;
mask_width_A0_2=100;
mask_width_S0_1=50;
mask_width_S0_2=80;
mask_width_SH0_1=30;
mask_width_SH0_2=50;
% Gauss window paramiters for blurring mask
Gauss_win_size = 10;
Gauss_sigma = 2; 
fprintf('Computing 3D dispersion relation: %s\n', modelname);
%% load optimized constants
output_name = [model_input_path,filesep,num2str(test_case),'output'];
load(output_name); % 'C11','C12','C13','C22','C23','C33','C44','C55','C66','rho','ObjVal'
%% dispersion curves - faster algorithm - wavenumber sweep
% [wavenumber1,CG1,FREQ1] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
% save([dispersion_output_path,'dispersion3D_uni1'],'wavenumber1','CG1','FREQ1');
disp('loading dispersion curves');
load([dispersion_output_path,'dispersion3D_uni1']);% 'wavenumber1','CG1','FREQ1' - dispersion curves for optimized constants
clear CG1;
% sanity check - plot A0 mode
% figure;
% for k=1:length(beta)   
%     plot(FREQ1(1,:,k),wavenumber1(:,k))
%     hold on;
% end
%%
% for other modes than A0, root sorting is needed
%%

%%
disp('loading wavefield data');
wavefield_name='333x333p_100kHz_5HC_10Vpp_x10_pzt';
load(['/pkudela_odroid_laser/aidd/data/raw/exp/L1_S2_B/',wavefield_name]);% Data, WL, time
Data=rot90(Data);
[m,n,nft]=size(Data);
Width=WL(1);
Length=WL(2);
disp('Transform to wavenumber-wavenumber-frequency domain');
[KXKYF,kx_vec,ky_vec,f_vec] = spatial_to_wavenumber_wavefield_full(Data,Length,Width,time);
% interpolate data at selected frequencies f_vec
disp('wavenumber interpolation');
wavenumber1q_A0=zeros(length(f_vec),length(beta));
wavenumber1q_S0=zeros(length(f_vec),length(beta));
wavenumber1q_SH0=zeros(length(f_vec),length(beta));
number_of_modes = 7;
for k=1:length(beta)   
    [k,length(beta)]
    [FREQ_new] = mode_tracing_new(squeeze(FREQ1(:,:,k)),wavenumber1(:,k),number_of_modes);
    wavenumber1q_A0(:,k) = interp1(FREQ_new(:,1)',wavenumber1(:,k),f_vec,'spline'); % A0
    wavenumber1q_S0(:,k) = interp1(FREQ_new(:,3)',wavenumber1(:,k),f_vec,'spline'); % S0
    wavenumber1q_SH0(:,k) = interp1(FREQ_new(:,2)',wavenumber1(:,k),f_vec,'spline'); % SH0
end
% correction of numerical errors
wavenumber1q_A0(1,:)=abs(wavenumber1q_A0(1,:));
wavenumber1q_S0(1,:)=abs(wavenumber1q_S0(1,:));
wavenumber1q_SH0(1,:)=abs(wavenumber1q_SH0(1,:));
temp = wavenumber1q_SH0(296:512,1:2);
wavenumber1q_SH0(296:512,1:2)=wavenumber1q_S0(296:512,1:2);
wavenumber1q_S0(296:512,1:2) = temp;

%for k=1:512 polarplot(beta*pi/180,wavenumber1q_SH0(k,:),'g');hold on;polarplot(beta*pi/180,wavenumber1q_S0(k,:),'b');polarplot(beta*pi/180,wavenumber1q_A0(k,:),'r');pause;clf;end;

%% create mask
disp('compute mask');

% [rad/m] mask width lineary varies with frequencies (actually this is half width)
% mask at lower frequencies is narrower
% mask_width = [linspace(0,mask_width_1,5),linspace(mask_width_1,mask_width_2,length(f_vec)-5)];

% mask A0 (Cartesian coordinates) - loop form (slow)
% mask_A0=zeros(length(kx_vec),length(ky_vec),length(f_vec));
% for fn=1:length(f_vec) % loop over frequencies (mask for each frequency bin)
%     [fn,length(f_vec)]
%     w_l = squeeze(wavenumber1q_A0(fn, :)-mask_width(fn));
%     w_u = squeeze(wavenumber1q_A0(fn, :)+mask_width(fn));
%     for j=1:length(ky_vec)
%         for i=1:length(kx_vec)
%             k=sqrt(kx_vec(i)^2+ky_vec(j)^2);
%             if(k==0) 
%                 b=0;
%             else
%                 if(kx_vec(i)>ky_vec(j))
%                     b=acos(kx_vec(i)/k)*180/pi;
%                 else
%                     b=asin(ky_vec(j)/k)*180/pi;
%                 end
%             end
%             [A,I]=min(abs(beta-b));
%             if(k>w_l(I) && k<w_u(I))
%                 mask_A0(j,i,fn)=1;
%             end
%         end
%     end
% end
%% vectorized form of mask calculation (fast)
[kx_grid,ky_grid]=meshgrid(kx_vec,ky_vec);
k_grid=sqrt(kx_grid.^2+ky_grid.^2);% wavenumber values on interpolated grid
b_grid= zeros(length(kx_vec),length(ky_vec));% angle values on interpolated grid
I1=(kx_grid>ky_grid);
b_grid(I1)=acos(kx_grid(I1)./k_grid(I1))*180/pi;
I2=(kx_grid<=ky_grid);I2(1,1)=0;
b_grid(I2)=asin(ky_grid(I2)./k_grid(I2))*180/pi;
k_gridp=repmat(k_grid,1,1,length(f_vec));
ind=zeros(length(kx_vec),length(ky_vec));
for j=1:length(ky_vec)
    for i=1:length(kx_vec)
        [A,I]=min(abs(beta-b_grid(i,j)));
        ind(i,j)=I;
    end
end
%% mask A0 
    mask_width = [linspace(0,mask_width_A0_1,5),linspace(mask_width_A0_1,mask_width_A0_2,length(f_vec)-5)];
    wavenumber_lower_bound = (wavenumber1q_A0 - repmat(mask_width',1,length(beta)))'; %[length(beta),length(f_vec)]
    wavenumber_upper_bound = (wavenumber1q_A0 + repmat(mask_width',1,length(beta)))';
   
    J1=(k_gridp > reshape(wavenumber_lower_bound(reshape(ind,[],1),:),length(kx_vec),length(ky_vec),length(f_vec)));
    J2=(k_gridp < reshape(wavenumber_upper_bound(reshape(ind,[],1),:),length(kx_vec),length(ky_vec),length(f_vec)));
    mask_A0 = J1.*J2;
    %fn=85;figure;surf(mask_A0(:,:,fn));shading interp; view(2);axis equal;
    
    % mask for all quadrants (mirroring)
    mask_A0_4quadrants=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
    mask_A0_4quadrants(length(kx_vec)+1:2*length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=mask_A0;
    mask_A0_4quadrants(1:length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=flipud(mask_A0);
    mask_A0_4quadrants(length(kx_vec)+1:2*length(kx_vec),1:length(ky_vec),:)=fliplr(mask_A0);
    mask_A0_4quadrants(1:length(kx_vec),1:length(ky_vec),:)=rot90(mask_A0,2);
    % blur mask
    [H] = Gauss(Gauss_win_size,Gauss_sigma);  
    FilterMaskBlur=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
    for k=1:length(f_vec)
        FilterMaskBlur(:,:,k) = filter2(H,mask_A0_4quadrants(:,:,k));
    end
    % mask for all quadrants symmetric in frequencies
    mask_A0_4quadrants_sym=zeros(2*length(kx_vec),2*length(ky_vec),2*length(f_vec));
    mask_A0_4quadrants_sym(:,:,1:length(f_vec))=flip(FilterMaskBlur,3);
    mask_A0_4quadrants_sym(:,:,length(f_vec)+1:2*length(f_vec))=FilterMaskBlur;
    clear mask_A0_4quadrants mask_A0; 
%% mask S0
    mask_width = [linspace(0,mask_width_S0_1,5),linspace(mask_width_S0_1,mask_width_S0_2,length(f_vec)-5)];
    wavenumber_lower_bound = (wavenumber1q_S0 - repmat(mask_width',1,length(beta)))'; %[length(beta),length(f_vec)]
    wavenumber_upper_bound = (wavenumber1q_S0 + repmat(mask_width',1,length(beta)))';
    J1=(k_gridp > reshape(wavenumber_lower_bound(reshape(ind,[],1),:),length(kx_vec),length(ky_vec),length(f_vec)));
    J2=(k_gridp < reshape(wavenumber_upper_bound(reshape(ind,[],1),:),length(kx_vec),length(ky_vec),length(f_vec)));
    mask_S0 = J1.*J2;
    %figure;surf(mask_A0(:,:,fn));shading interp; view(2);axis equal;
    
    % mask for all quadrants (mirroring)
    mask_S0_4quadrants=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
    mask_S0_4quadrants(length(kx_vec)+1:2*length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=mask_S0;
    mask_S0_4quadrants(1:length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=flipud(mask_S0);
    mask_S0_4quadrants(length(kx_vec)+1:2*length(kx_vec),1:length(ky_vec),:)=fliplr(mask_S0);
    mask_S0_4quadrants(1:length(kx_vec),1:length(ky_vec),:)=rot90(mask_S0,2);
    % blur mask
    [H] = Gauss(Gauss_win_size,Gauss_sigma);  
    FilterMaskBlur=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
    for k=1:length(f_vec)
        FilterMaskBlur(:,:,k) = filter2(H,mask_S0_4quadrants(:,:,k));
    end
    % mask for all quadrants symmetric in frequencies
    mask_S0_4quadrants_sym=zeros(2*length(kx_vec),2*length(ky_vec),2*length(f_vec));
    mask_S0_4quadrants_sym(:,:,1:length(f_vec))=flip(FilterMaskBlur,3);
    mask_S0_4quadrants_sym(:,:,length(f_vec)+1:2*length(f_vec))=FilterMaskBlur;
    clear mask_S0_4quadrants; 
%% mask SH0
    mask_width = [linspace(0,mask_width_SH0_1,5),linspace(mask_width_SH0_1,mask_width_SH0_2,length(f_vec)-5)];
    wavenumber_lower_bound = (wavenumber1q_SH0 - repmat(mask_width',1,length(beta)))'; %[length(beta),length(f_vec)]
    wavenumber_upper_bound = (wavenumber1q_SH0 + repmat(mask_width',1,length(beta)))';
    J1=(k_gridp > reshape(wavenumber_lower_bound(reshape(ind,[],1),:),length(kx_vec),length(ky_vec),length(f_vec)));
    J2=(k_gridp < reshape(wavenumber_upper_bound(reshape(ind,[],1),:),length(kx_vec),length(ky_vec),length(f_vec)));
    mask_SH0 = J1.*J2;
    
    % exclude common mask elements with S0
    mask_S0_SH0_intersect = mask_S0 .* mask_SH0;
    mask_SH0 = mask_SH0 - mask_S0_SH0_intersect;
    %figure;surf(mask_SH0(:,:,fn));shading interp; view(2);axis equal;
    
    % mask for all quadrants (mirroring)
    mask_SH0_4quadrants=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
    mask_SH0_4quadrants(length(kx_vec)+1:2*length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=mask_SH0;
    mask_SH0_4quadrants(1:length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=flipud(mask_SH0);
    mask_SH0_4quadrants(length(kx_vec)+1:2*length(kx_vec),1:length(ky_vec),:)=fliplr(mask_SH0);
    mask_SH0_4quadrants(1:length(kx_vec),1:length(ky_vec),:)=rot90(mask_SH0,2);
    % blur mask
    [H] = Gauss(Gauss_win_size,Gauss_sigma);  
    FilterMaskBlur=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
    for k=1:length(f_vec)
        FilterMaskBlur(:,:,k) = filter2(H,mask_SH0_4quadrants(:,:,k));
    end
    % mask for all quadrants symmetric in frequencies
    mask_SH0_4quadrants_sym=zeros(2*length(kx_vec),2*length(ky_vec),2*length(f_vec));
    mask_SH0_4quadrants_sym(:,:,1:length(f_vec))=flip(FilterMaskBlur,3);
    mask_SH0_4quadrants_sym(:,:,length(f_vec)+1:2*length(f_vec))=FilterMaskBlur;   
    
    clear mask_SH0_4quadrants mask_SH0 mask_S0; 
% clear unnecessary matrices
clear k_grid k_gridp kx_grid ky_grid;
%% apply mask
disp('processing and inverse Fourier transform');
KXKYF_A0 = KXKYF.*mask_A0_4quadrants_sym;
KXKYF_S0 = KXKYF.*mask_S0_4quadrants_sym;
KXKYF_SH0 = KXKYF.*mask_SH0_4quadrants_sym;
%% inverse Fourier transform for pure A0, S0, SH0 wavefield
W_A0 = ifftn(ifftshift(KXKYF_A0),'symmetric'); % wavefield A0
W_S0 = ifftn(ifftshift(KXKYF_S0),'symmetric'); % wavefield S0 
W_SH0 = ifftn(ifftshift(KXKYF_SH0),'symmetric'); % wavefield SH0
%% plotting
disp('plotting');
fn=50; % f_vec(80)=9.8943e+04; %[Hz]
figure;surf(mask_A0_4quadrants_sym(:,:,length(f_vec)+fn));shading interp;view(2);axis square;xlim([1 2*length(kx_vec)]);ylim([1 2*length(ky_vec)]);
figure;surf((squeeze(mask_A0_4quadrants_sym(end/2+1:end,end/2,length(f_vec)+1:end))));shading interp;view(2);xlim([1 length(kx_vec)]);ylim([1 length(f_vec)]);

figure;surf(mask_S0_4quadrants_sym(:,:,length(f_vec)+fn));shading interp;view(2);axis square;xlim([1 2*length(kx_vec)]);ylim([1 2*length(ky_vec)]);
figure;surf((squeeze(mask_S0_4quadrants_sym(end/2+1:end,end/2,length(f_vec)+1:end))));shading interp;view(2);xlim([1 length(kx_vec)]);ylim([1 length(f_vec)]);

figure;surf(mask_SH0_4quadrants_sym(:,:,length(f_vec)+fn));shading interp;view(2);axis square;xlim([1 2*length(kx_vec)]);ylim([1 2*length(ky_vec)]);
figure;surf((squeeze(mask_SH0_4quadrants_sym(end/2+1:end,end/2,length(f_vec)+1:end))));shading interp;view(2);xlim([1 length(kx_vec)]);ylim([1 length(f_vec)]);


% figure;surf(mask_S0_SH0_intersect(:,:,length(f_vec)+fn));shading interp;view(2);axis square;xlim([1 2*length(kx_vec)]);ylim([1 2*length(ky_vec)]);
% figure;surf(mask_SH0_4quadrants_sym(:,:,length(f_vec)+fn));shading interp;view(2);axis square;xlim([1 2*length(kx_vec)]);ylim([1 2*length(ky_vec)]);

%figure;surf(squeeze(abs(KXKYF(2:end,2:end,512+fn))));shading interp;view(2);axis equal;xlim([1 2*length(kx_vec)]); ylim([1 2*length(ky_vec)]);axis off;
%figure;surf(squeeze(abs(KXKYF(2:end,2:end,512-fn))));shading interp;view(2);axis equal;axis equal;xlim([1 2*length(kx_vec)]); ylim([1 2*length(ky_vec)]);axis off;
%figure;surf(squeeze(mask_A0_4quadrants_sym(2:end,2:end,512+fn)));shading interp;view(2);axis equal;axis equal;xlim([1 2*length(kx_vec)]); ylim([1 2*length(ky_vec)]);axis off;
%figure;surf(squeeze(abs(KXKYF_A0(2:end,2:end,512+fn))));shading interp;view(2);axis equal;axis equal;xlim([1 2*length(kx_vec)]); ylim([1 2*length(ky_vec)]);axis off;
%% plot figures - checking results
filename_A0 = [wavefield_name,'_A0'];
filename_S0 = [wavefield_name,'_S0'];
filename_SH0 = [wavefield_name,'_SH0'];
for fn=[50,100,150,200,250,300]
    fh=figure('Menu','none','ToolBar','none');
        surf(squeeze((real(W_A0(1:m,1:n,fn))))); shading interp;view(2);colormap jet;xlim([1 m]); ylim([1 n]);axis off;
        figfilename = [filename_A0,'_frame_',num2str(fn)];
        set(gcf,'Color','w');set(gca,'FontName','Times');set(gca,'Fontsize',10,'linewidth',1); set(gcf,'Renderer','zbuffer');
        fig = gcf;set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
        % remove unnecessary white space
        set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
        set(gcf,'PaperPositionMode','auto');
        print([figure_output_path,figfilename],'-dpng', '-r600'); 
    figure('Menu','none','ToolBar','none');surf(squeeze((real(W_S0(1:m,1:n,fn)))));shading interp;view(2);colormap jet;xlim([1 m]); ylim([1 n]);axis off;
        figfilename = [filename_S0,'_frame_',num2str(fn)];
        set(gcf,'Color','w');set(gca,'FontName','Times');set(gca,'Fontsize',10,'linewidth',1); set(gcf,'Renderer','zbuffer');
        fig = gcf;set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
        % remove unnecessary white space
        set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
        set(gcf,'PaperPositionMode','auto');
        print([figure_output_path,figfilename],'-dpng', '-r600'); 
    figure('Menu','none','ToolBar','none');surf(squeeze((real(W_SH0(1:m,1:n,fn)))));shading interp;view(2);colormap jet;xlim([1 m]); ylim([1 n]);axis off;
        figfilename = [filename_SH0,'_frame_',num2str(fn)];
        set(gcf,'Color','w');set(gca,'FontName','Times');set(gca,'Fontsize',10,'linewidth',1); set(gcf,'Renderer','zbuffer');
        fig = gcf;set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
        % remove unnecessary white space
        set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
        set(gcf,'PaperPositionMode','auto');
        print([figure_output_path,figfilename],'-dpng', '-r600'); 
    figure('Menu','none','ToolBar','none');surf(squeeze(((Data(:,:,fn)))));shading interp;view(2);colormap jet;xlim([1 m]); ylim([1 n]);axis off;
        figfilename = [wavefield_name,'_frame_',num2str(fn)];
        set(gcf,'Color','w');set(gca,'FontName','Times');set(gca,'Fontsize',10,'linewidth',1); set(gcf,'Renderer','zbuffer');
        fig = gcf;set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
        % remove unnecessary white space
        set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
        set(gcf,'PaperPositionMode','auto');
        print([figure_output_path,figfilename],'-dpng', '-r600'); 
end

%% saving data
Data_A0 = real(W_A0(1:m,1:n,:));
Data_S0 = real(W_S0(1:m,1:n,:));
Data_SH0 = real(W_SH0(1:m,1:n,:));
save([output_path,filename_A0],'Data_A0','WL','time');
save([output_path,filename_S0],'Data_S0','WL','time');
save([output_path,filename_SH0],'Data_SH0','WL','time');
