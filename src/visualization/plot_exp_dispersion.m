% plot experimental dispersion relations depending on angle

clear all; close all;
% load projectroot path
load project_paths projectroot src_path;

% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
% and conversion to polar coordinate system
filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_'; 

% load experimental data file
disp('loading data ...');
load([data_path,filename]); % Data_polar x y wavenumber_max fmax

[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);
for j=1:number_of_angles
    figure(j);
%     surf(squeeze(abs(Data_polar(j,2:end,2:end))));
%     %surf(squeeze(abs(real(Data_polar(j,2:end,2:end)))));
%     surf(squeeze(abs(imag(Data_polar(j,2:end,2:end)))));
%     shading interp; view(2);colormap jet
    kvec = linspace(0,wavenumber_max(j),number_of_wavenumber_points);
    set(gcf,'Color','w');
    imagesc(fvec/1e3, kvec, squeeze(abs(Data_polar(j,2:end,2:end))));   set(gca,'YDir','normal'); 
    axis([0 600 0 3000]);
    set(gca,'Fontsize',14)
    xlabel('f [kHz]')
    ylabel('k [1/m]')
    colormap jet; 
%axis tight; 
    caxis([0 max(caxis)/3]);
    
end


j=1; % beta = 0
