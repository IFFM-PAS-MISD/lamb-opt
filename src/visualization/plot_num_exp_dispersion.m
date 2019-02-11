% plot numerical dispersion curves on top of experimental data

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;

% create path to the numerical model data folder
foldername = 'SASE';
modelname = 'SASE1';

% create output path
output_path = prepare_figure_paths(foldername,modelname);

% create path to the numerical raw data folder
model_input_path = fullfile( projectroot, 'data','raw','num',foldername,[modelname,'_out'], filesep );
% create path to the experimental processed data folder
exp_input_data_path=fullfile( projectroot, 'data','processed','exp', filesep );

% filename of data to be processed
% full field measurements after 3D FFT transform (1st quarter)
% and conversion to polar coordinate system
filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_'; 

% load experimental data file
disp('loading data ...');
load([exp_input_data_path,filename]); % Data_polar x y wavenumber_max fmax
for test_case=[25,29,121]
% load numerical data file
%test_case = 121; % worst objective function score
%test_case = 25; % best objective function score for all modes
%test_case = 29; % best objective function score for A0 mode only

input_name = [model_input_path,num2str(test_case),'output'];
load(input_name); % FREQ CG wavenumber

%% START PLOTTING
[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);
beta=0:90/(number_of_angles-1):90;
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
    fvec1=squeeze(FREQ(1,:,j)); % mode 1, angle j
    fvec2=squeeze(FREQ(2,:,j)); % mode 2, angle j
    fvec3=squeeze(FREQ(3,:,j)); % mode 3, angle j
    fvec4=squeeze(FREQ(4,:,j)); % mode 4, angle j
    kvec=squeeze(wavenumber(:,j)); % angle j
    plot(fvec1(2:end)/1e3,kvec(2:end),'y');
    plot(fvec2(2:end)/1e3,kvec(2:end),'y');
    plot(fvec3(2:end)/1e3,kvec(2:end),'y');
    plot(fvec4(2:end)/1e3,kvec(2:end),'y');
    title(['angle ', num2str(beta(j)),' deg']);
    figfilename = [num2str(test_case),'_angle_',num2str(beta(j)),'_num_exp_dispersion'];
    print('-dpng', '-r300',[output_path,figfilename]); 
end
close all;
end
%% END PLOTTING

