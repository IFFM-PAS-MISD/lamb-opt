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

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
input_file = 1; % experimental data
% load experimental data file
exp_filename = {'interim_499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni_KXKYF',... % 1 small area unidirectional
                            'interim_499x499p_chp200_x40_6Vpp_250Hz_uni_KXKYF'};         % 2 large area unidirectional
load([exp_input_data_path,exp_filename{input_file}]); % 'KXKYF_','f_vec','kx_vec', 'ky_vec' 
% filename of parameter data
filename = {'polar_f_interim_499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni_KXKYF_param',...
                    'polar_f_interim_499x499p_chp200_x40_6Vpp_250Hz_uni_KXKYF_param'}; 
load([exp_param_data_path,filename{input_file}]); % beta number_of_wavenumber_points selected_frequency_index selected_frequencies ...
fmin = f_vec(2); % first frequency is 0
fmax = f_vec(end);
number_of_frequency_points = length(f_vec)-1;
%% Input for SASE
ht = 2.85/1000; % [m] laminate total thickness; unidirectional
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
wavenumber_max = zeros(length(beta),1)+3.43e3; % maximal wavenumbers [1/m]
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

Data = KXKYF_(:,:,selected_frequency_index);
clear KXKYF_;
[number_of_wavenumber_points_x,number_of_wavenumber_points_y,~] = size(Data);



% if(~radians)
%    wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
% end
fprintf('Computing 3D dispersion relation: %s\n', modelname);
%% load optimized constants
output_name = [model_input_path,filesep,num2str(test_case),'output'];
load(output_name); % 'C11','C12','C13','C22','C23','C33','C44','C55','C66','rho','ObjVal'

%% compute dispersion curves
%[wavenumber,CG,FREQ] = main_SASE2(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,fmin,fmax,number_of_frequency_points,beta,stack_dir,np,nele_layer);
%save([output_path,'dispersion3D_uni2'],'wavenumber','CG','FREQ');
%load([output_path,'dispersion3D_uni2']);% 'wavenumber','CG','FREQ'
%% sort dispersion curve roots
% number_of_modes=7;
% wavenumber_sorted=zeros(number_of_modes,number_of_frequency_points,length(beta));
% % loop over angles
% for k=1:length(beta)
%     [k,length(beta)]
%     [wavenumber_new] = mode_tracing_new2(FREQ(:,k),squeeze(wavenumber(:,:,k)),number_of_modes);
%     wavenumber_sorted(:,:,k) = wavenumber_new';
% end
% % find A0 mode
% selected_freq_no = 250;expected_wavenumber=1500;
% wavenumber_A0=zeros(number_of_frequency_points,length(beta));
% for k=1:length(beta)
%     [A,I]=min(abs(wavenumber_sorted(:,selected_freq_no,k)-expected_wavenumber));
%     wavenumber_A0(:,k) = wavenumber_sorted(I,:,k);
% end
% % sanity check - plot
% figure;
% for k=1:length(beta)   
%     [k]
%     plot(wavenumber_A0(:,k));
%     %hold on;
%     pause;clf;
% end
%% dispersion curves - faster algorithm - wavenumber sweep
% [wavenumber1,CG1,FREQ1] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
% save([output_path,'dispersion3D_uni1'],'wavenumber1','CG1','FREQ1');
load([output_path,'dispersion3D_uni1']);% 'wavenumber1','CG1','FREQ1'
% sanity check - plot A0 mode
% figure;
% for k=1:length(beta)   
%     plot(FREQ1(1,:,k),wavenumber1(:,k))
%     hold on;
% end
%%
% interpolate data at selected frequencies f_vec
wavenumber1q_A0=zeros(length(f_vec),length(beta));
for k=1:length(beta)   
    wavenumber1q_A0(:,k) = interp1(FREQ1(1,:,k),wavenumber1(:,k),f_vec,'spline'); % for A0 mode only
end


 % transform numerical dispersion surface from (k-beta) polar coordinate system to (kx-ky) Cartesian coordinate system
mask_width = linspace(20,80,length(f_vec)); % [rad/m] mask width lineary varies with frequencies (actually this is half width)
mask_width(1:10)=linspace(0,10,10); % mask width at lower frequencies is narrower
% wx_lower = zeros(length(f_vec),length(beta)); 
% wy_lower = zeros(length(f_vec),length(beta));
% wx_center = zeros(length(f_vec),length(beta));
% wy_center = zeros(length(f_vec),length(beta));
% wx_upper = zeros(length(f_vec),length(beta));
% wy_upper = zeros(length(f_vec),length(beta));
% for j=2:length(f_vec)
%     wx_lower(j,:) = squeeze(wavenumber1q_A0(j, :)-mask_width(j)).*cos(beta*pi/180); % % lower bound of the wavenumber of the mask
%     wy_lower(j,:) = squeeze(wavenumber1q_A0(j, :)-mask_width(j)).*sin(beta*pi/180); % % lower bound of the wavenumber of the mask
%     wx_center(j,:) = squeeze(wavenumber1q_A0(j, :)).*cos(beta*pi/180); % 
%     wy_center(j,:) = squeeze(wavenumber1q_A0(j, :)).*sin(beta*pi/180); % 
%     wx_upper(j,:) = squeeze(wavenumber1q_A0(j, :)+mask_width(j)).*cos(beta*pi/180); % % upper bound of the wavenumber of the mask
%     wy_upper(j,:) = squeeze(wavenumber1q_A0(j, :)+mask_width(j)).*sin(beta*pi/180); % % upper bound of the wavenumber of the mask
% end
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

% mirroring just for checking
% KXKYF=zeros(2*length(kx_vec),2*length(ky_vec),length(f_vec));
% KXKYF(length(kx_vec)+1:2*length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=KXKYF_A0;
% KXKYF(1:length(kx_vec),length(ky_vec)+1:2*length(ky_vec),:)=flipud(KXKYF_A0);
% KXKYF(length(kx_vec)+1:2*length(kx_vec),1:length(ky_vec),:)=fliplr(KXKYF_A0);
% KXKYF(1:length(kx_vec),1:length(ky_vec),:)=rot90(KXKYF_A0,2);
% 
% KXKYF2=flip(KXKYF,3);
% figure;surf(squeeze(abs(KXKYF(1:end,1:end,fn))));shading interp;view(2);axis equal;
% 
% 
% S=zeros(2*length(kx_vec),2*length(ky_vec),2*length(f_vec));
% S(:,:,1:length(f_vec))=KXKYF2;
% S(:,:,length(f_vec)+1:2*length(f_vec))=KXKYF;
% W = ifftn(ifftshift(S,3),'symmetric');
% figure;surf(W(:,:,200));shading interp;view(2);axis equal;
% W = ifftn(ifftshift(S));
% figure;surf(real(W(:,:,200)));shading interp;view(2);axis equal;
% W = ifftn(ifftshift(KXKYF));
% 
% W = ifftn(ifftshift(KXKYF_),'symmetric');
return;


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