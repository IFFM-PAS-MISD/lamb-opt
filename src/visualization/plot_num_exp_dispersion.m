% plot numerical dispersion curves on top of experimental data

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;

% create path to the numerical model data folder
foldername = 'SASE';
modelname = 'SASE1';
data_process_type = 'raw';
data_origin = 'num';
model_output_path = prepare_model_paths(data_process_type,data_origin,foldername,modelname);

% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
% and conversion to polar coordinate system
filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_'; 

% load experimental data file
disp('loading data ...');
load([data_path,filename]); % Data_polar x y wavenumber_max fmax
% load numerical data file
test_case = 121; % input
output_name = [model_output_path,filesep,'output',num2str(test_case)];
load(output_name); % FREQ CG wavenumber

%% START PLOTTING
[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);

for j=1:number_of_angles % beta

figure(j)
set(gcf,'Color','w');
%imagesc(fvec(2:end)/1e3, wavenumber(2:end,j)/2, squeeze(abs(Data_polar(j,2:end,2:end)))); 
imagesc(fvec(2:end)/1e3, wavenumber(2:end,j), squeeze(abs(Data_polar(j,2:end,2:end)))); 
set(gca,'YDir','normal'); 
%axis([0 600 0 2000]);
%axis([0 350 0 2000]);
axis([0 350 0 min(wavenumber_max)]);
set(gca,'Fontsize',14)
xlabel('f [kHz]')
ylabel('k [rd/m]')
colormap jet; 
%axis tight; 
caxis([0 max(caxis)/3]); 

hold on;
%plot(freq(1:25000)/1e3,wavenumber(1:25000),'y');
fvec1=squeeze(FREQ(1,:,j)); % mode 1, angle j
fvec2=squeeze(FREQ(2,:,j)); % mode 1, angle j
fvec3=squeeze(FREQ(3,:,j)); % mode 1, angle j
fvec4=squeeze(FREQ(4,:,j)); % mode 1, angle j
kvec=squeeze(wavenumber(:,j)); % angle j
plot(fvec1(2:end)/1e3,kvec(2:end),'y');
plot(fvec2(2:end)/1e3,kvec(2:end),'y');
plot(fvec3(2:end)/1e3,kvec(2:end),'y');
plot(fvec4(2:end)/1e3,kvec(2:end),'y');

end

%% END PLOTTING

