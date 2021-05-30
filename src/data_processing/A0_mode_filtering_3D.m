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
test_case=2; % numerical data
% create output path
output_path = prepare_figure_paths(modelfolder,modelname);

% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',modelfolder,[modelname,'_out'], filesep );
% create path to the experimental interim data folder
exp_input_data_path=fullfile( projectroot, 'data','interim','exp', filesep );
% create path to the experimental processed data folder
exp_param_data_path=fullfile( projectroot, 'data','processed','exp', filesep );
% create path to the experimental wavefield data folder
exp_wavefield_data_path=fullfile( projectroot, 'data','raw','exp', filesep );
exp_wavefield_filename='499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni';
load([exp_wavefield_data_path,exp_wavefield_filename]);

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
input_file = 1; % experimental data
% load experimental data file
exp_filename = {'interim_499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni_KXKYF',... % 1 small area unidirectional
                            'interim_499x499p_chp200_x40_6Vpp_250Hz_uni_KXKYF'};         % 2 large area unidirectional
load([exp_input_data_path,exp_filename{input_file}]); % 'KXKYF_','f_vec','kx_vec', 'ky_vec' 
% filename of parameter data
% filename = {'polar_f_interim_499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni_KXKYF_param',...
%                     'polar_f_interim_499x499p_chp200_x40_6Vpp_250Hz_uni_KXKYF_param'}; 
% load([exp_param_data_path,filename{input_file}]); % beta number_of_wavenumber_points selected_frequency_index selected_frequencies ...

%% Input for SASE
ht = 2.85/1000; % [m] laminate total thickness; unidirectional
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
no_of_angles=512;
beta = 0:90/(no_of_angles-1):90;
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
wavenumber_max = zeros(length(beta),1)+3.43e3; % maximal wavenumbers [1/m]
number_of_wavenumber_points=length(kx_vec);
layup = [90 90 90 90 90 90 90 90];

nlayers = length(layup);

h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers; new plain weave
% Stacking direction
stack_dir = 1;
%%
% known parameters
m=6.46; % total mass of the specimen [kg]
V=1.2*1.2*ht; % specimen volume [m^3]
rho = m/V;
%%
%surf(squeeze(abs(KXKYF_(2:end,20,2:end))));shading interp; view(2);axis square;colormap jet;


% if(~radians)
%    wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
% end
fprintf('Computing 3D dispersion relation: %s\n', modelname);
%% load optimized constants
output_name = [model_input_path,filesep,num2str(test_case),'output'];
load(output_name); % 'C11','C12','C13','C22','C23','C33','C44','C55','C66','rho','ObjVal'
%% dispersion curves - faster algorithm - wavenumber sweep
% [wavenumber1,CG1,FREQ1] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
% save([output_path,'dispersion3D_uni1'],'wavenumber1','CG1','FREQ1');
load([output_path,'dispersion3D_uni1']);% 'wavenumber1','CG1','FREQ1' - dispersion curves for optimized constants
% sanity check - plot A0 mode
% figure;
% for k=1:length(beta)   
%     plot(FREQ1(1,:,k),wavenumber1(:,k))
%     hold on;
% end
%%
% for other modes than A0, root sorting is needed
%%
% interpolate data at selected frequencies f_vec
wavenumber1q_A0=zeros(length(f_vec),length(beta));
for k=1:length(beta)   
    wavenumber1q_A0(:,k) = interp1(FREQ1(1,:,k),wavenumber1(:,k),f_vec,'spline'); % for A0 mode only
end

%% create mask
mask_width = linspace(20,80,length(f_vec)); % [rad/m] mask width lineary varies with frequencies (actually this is half width)
mask_width(1:10)=linspace(0,10,10); % mask width at lower frequencies is narrower

% mask A0 (Cartesian coordinates) - needs to be done in vectorized form
mask_A0=zeros(length(kx_vec),length(ky_vec),length(f_vec));
for fn=1:length(f_vec) % loop over frequencies (mask for each frequency bin)
    [fn,length(f_vec)]
    w_l = squeeze(wavenumber1q_A0(fn, :)-mask_width(fn));
    w_u = squeeze(wavenumber1q_A0(fn, :)+mask_width(fn));
    for j=1:length(ky_vec)
        for i=1:length(kx_vec)
            k=sqrt(kx_vec(i)^2+ky_vec(j)^2);
            if(k==0) 
                b=0;
            else
                if(kx_vec(i)>ky_vec(j))
                    b=acos(kx_vec(i)/k)*180/pi;
                else
                    b=asin(ky_vec(j)/k)*180/pi;
                end
            end
            [A,I]=min(abs(beta-b));
            if(k>w_l(I) && k<w_u(I))
                mask_A0(j,i,fn)=1;
            end
        end
    end
end

% sanity check
fn=100;
figure;surf(squeeze(abs(KXKYF_(2:end,2:end,fn))));shading interp;view(2);axis equal;
figure;surf(mask_A0(2:end,2:end,fn));shading interp;view(2);axis equal;
KXKYF_A0 = KXKYF_.*mask_A0;
figure;surf(squeeze(abs(KXKYF_A0(2:end,2:end,fn))));shading interp;view(2);axis equal;
% mask for all quadrants (mirroring)
mask_A0_4quadrants=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
mask_A0_4quadrants(length(kx_vec)+1:2*length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=mask_A0;
mask_A0_4quadrants(1:length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=flipud(mask_A0);
mask_A0_4quadrants(length(kx_vec)+1:2*length(kx_vec),1:length(ky_vec),:)=fliplr(mask_A0);
mask_A0_4quadrants(1:length(kx_vec),1:length(ky_vec),:)=rot90(mask_A0,2);
figure;surf(mask_A0_4quadrants(1:end,1:end,fn));shading interp;view(2);axis equal;
%% blur mask
% to be implemented
% mirroring just for checking
KXKYF=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
KXKYF(length(kx_vec)+1:2*length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=KXKYF_A0;
KXKYF(1:length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=flipud(KXKYF_A0);
KXKYF(length(kx_vec)+1:2*length(kx_vec),1:length(ky_vec),:)=fliplr(KXKYF_A0);
KXKYF(1:length(kx_vec),1:length(ky_vec),:)=rot90(KXKYF_A0,2);
figure;surf(squeeze(abs(KXKYF(2:end,2:end,fn))));shading interp;view(2);axis equal;
% 

% 
S=zeros(2*length(kx_vec),2*length(ky_vec),2*length(f_vec));
S(:,:,1:length(f_vec))=flip(KXKYF,3);
S(:,:,length(f_vec)+1:2*length(f_vec))=KXKYF;
%W = ifftn(ifftshift(S)); % wavefield is on the first quarter of the marix (complex)
W = ifftn(ifftshift(S),'symmetric'); % wavefield is on the first quarter of the marix (always real)

figure;surf(squeeze((real(W(1:end/2,1:end/2,100)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,100)))));shading interp;view(2);axis equal;colormap jet;
return;
%%

load(['/pkudela_odroid_laser/aidd/data/raw/exp/L1_S2_B/333x333p_100kHz_5HC_10Vpp_x10_pzt']);% Data, WL, time
Data=rot90(Data);
[m,n,nft]=size(Data);
Width=WL(1);
Length=WL(2);
[KXKYF,kx_vec,ky_vec,f_vec] = spatial_to_wavenumber_wavefield_full(Data,Length,Width,time);
% interpolate data at selected frequencies f_vec
wavenumber1q_A0=zeros(length(f_vec),length(beta));
for k=1:length(beta)   
    wavenumber1q_A0(:,k) = interp1(FREQ1(1,:,k),wavenumber1(:,k),f_vec,'spline'); % for A0 mode only
end
%% create mask
mask_width = linspace(20,80,length(f_vec)); % [rad/m] mask width lineary varies with frequencies (actually this is half width)
mask_width(1:10)=linspace(0,10,10); % mask width at lower frequencies is narrower

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
% vectorized form of mask calculation (fast)
wavenumber_lower_bound = (wavenumber1q_A0 - repmat(mask_width',1,length(beta)))'; %[length(beta),length(f_vec)]
wavenumber_upper_bound = (wavenumber1q_A0 + repmat(mask_width',1,length(beta)))';
[kx_grid,ky_grid]=meshgrid(kx_vec,ky_vec);
k_grid=sqrt(kx_grid.^2+ky_grid.^2);% wavenumber values on interpolated grid
b_grid= zeros(length(kx_vec),length(ky_vec));% angle values on interpolated grid
I1=(kx_grid>ky_grid);
b_grid(I1)=acos(kx_grid(I1)./k_grid(I1))*180/pi;
I2=(kx_grid<=ky_grid);I2(1,1)=0;
b_grid(I2)=asin(ky_grid(I2)./k_grid(I2))*180/pi;
ind=zeros(length(kx_vec),length(ky_vec));
for j=1:length(ky_vec)
    for i=1:length(kx_vec)
        [A,I]=min(abs(beta-b_grid(i,j)));
        ind(i,j)=I;
    end
end
fn=85;
k_gridp=repmat(k_grid,1,1,length(f_vec));
J1=(k_gridp > reshape(wavenumber_lower_bound(reshape(ind,[],1),:),length(kx_vec),length(ky_vec),length(f_vec)));
J2=(k_gridp < reshape(wavenumber_upper_bound(reshape(ind,[],1),:),length(kx_vec),length(ky_vec),length(f_vec)));
mask_A0 = J1.*J2;
%figure;surf(mask_A0(:,:,fn));shading interp; view(2);axis equal;

%%
% mask for all quadrants (mirroring)
mask_A0_4quadrants=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
mask_A0_4quadrants(length(kx_vec)+1:2*length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=mask_A0;
mask_A0_4quadrants(1:length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=flipud(mask_A0);
mask_A0_4quadrants(length(kx_vec)+1:2*length(kx_vec),1:length(ky_vec),:)=fliplr(mask_A0);
mask_A0_4quadrants(1:length(kx_vec),1:length(ky_vec),:)=rot90(mask_A0,2);
% mask for all quadrants symmetric in frequencies
mask_A0_4quadrants_sym=zeros(2*length(kx_vec),2*length(ky_vec),2*length(f_vec));
mask_A0_4quadrants_sym(:,:,1:length(f_vec))=flip(mask_A0_4quadrants,3);
mask_A0_4quadrants_sym(:,:,length(f_vec)+1:2*length(f_vec))=mask_A0_4quadrants;
clear mask_A0_4quadrants;
%% blur mask
% to be implemented
figure;surf(mask_A0_4quadrants(1:end,1:end,fn));shading interp;view(2);axis equal;
% apply mask
KXKYF_A0 = KXKYF.*mask_A0_4quadrants_sym;
% check
fn=80; % f_vec(80)=9.8943e+04; %[Hz]
%figure;surf(squeeze(abs(KXKYF(2:end,2:end,512+fn))));shading interp;view(2);axis equal;
%figure;surf(squeeze(abs(KXKYF(2:end,2:end,512-fn))));shading interp;view(2);axis equal;
%figure;surf(squeeze(mask_A0_4quadrants_sym(2:end,2:end,512+fn)));shading interp;view(2);axis equal;
%figure;surf(squeeze(abs(KXKYF_A0(2:end,2:end,512+fn))));shading interp;view(2);axis equal;
% inverse Fourier transform for pure A0 wavefield
W_A0 = ifftn(ifftshift(KXKYF_A0),'symmetric'); % wavefield is on the first quarter of the marix (always real)
%% plot figures - checking results
%figure;surf(squeeze((real(W_A0(1:end/2,1:end/2,100)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze((real(W_A0(1:m,1:n,100)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,100)))));shading interp;view(2);axis equal;colormap jet;

figure;surf(squeeze((real(W_A0(1:m,1:n,150)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,150)))));shading interp;view(2);axis equal;colormap jet;

figure;surf(squeeze((real(W_A0(1:m,1:n,250)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,250)))));shading interp;view(2);axis equal;colormap jet;

figure;surf(squeeze((real(W_A0(1:m,1:n,350)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,350)))));shading interp;view(2);axis equal;colormap jet;

% inverse mask
mask_remaining_modes = ((-1)*mask_A0_4quadrants_sym)+1;
% apply mask
KXKYF_remaining_modes = KXKYF.*mask_remaining_modes;
% inverse Fourier transform for wavefield containing remaining nodes
W_remaining_modes = ifftn(ifftshift(KXKYF_remaining_modes),'symmetric'); % wavefield is on the first quarter of the marix (always real)
%% plot figures - checking results

figure;surf(squeeze((real(W_remaining_modes(1:m,1:n,100)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,100)))));shading interp;view(2);axis equal;colormap jet;

figure;surf(squeeze((real(W_remaining_modes(1:m,1:n,150)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,150)))));shading interp;view(2);axis equal;colormap jet;

figure;surf(squeeze((real(W_remaining_modes(1:m,1:n,250)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,250)))));shading interp;view(2);axis equal;colormap jet;

figure;surf(squeeze((real(W_remaining_modes(1:m,1:n,350)))));shading interp;view(2);axis equal;colormap jet;
figure;surf(squeeze(((Data(:,:,350)))));shading interp;view(2);axis equal;colormap jet;
return
%%
load project_paths projectroot src_path;
% input
overwrite = false; % allow overwriting existing results if true
% length and width of full wavefield area
Length = [0.455, 0.725];    
Width = [0.455, 0.725];
test_case=[1,2]; % select file numbers for processing

% create path to the experimental raw data folder
raw_data_path = fullfile( projectroot, 'data','raw','exp', filesep );

% create path to the experimental interim data folder
interim_data_path = fullfile( projectroot, 'data','interim','exp', filesep );

% filenames of data to be processed
% full field measurements
list = {'499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni', ...          % 1  Length = 0.455;Width = 0.455;           
        '499x499p_chp200_x40_6Vpp_250Hz_uni' };                % 2 Length = 0.724;Width = 0.727;
k=1;
filename = list{k};
load([raw_data_path,filename]); % Data Length Width time
[KXKYF,kx_vec,ky_vec,f_vec] = spatial_to_wavenumber_wavefield_full(Data,Length(k),Width(k),time);
%W = ifftn(ifftshift(KXKYF(:,:,:)),'symmetric');
W = ifftn(ifftshift(KXKYF(:,:,:)));
%W = ifftn(ifftshift(KXKYF(:,:,floor(end/2)+1:end)),'symmetric');% something is wrong
figure;surf(squeeze((W(1:end/2,1:end/2,100))));shading interp;view(2);axis equal;