% plot numerical dispersion curves on top of experimental data

clear all; close all;

%set(0,'defaulttextinterpreter','none');

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true

% figure parameters
% size 12cm by 8cm (1-column text)
%fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
fig_width = 7; fig_height = 5; 
% create path to the numerical model data folder
foldername = 'minigrant_raport_2';
radians = false;
% create output path
output_path = prepare_exp_figure_paths(foldername);
% length and width of full wavefield area
Length = 0.493;    
Width = 0.493;
test_case=[9:11,4:8,13,16:19,14,15,12,20,21,22]; % select file numbers for processing

% create path to the experimental raw data folder
raw_data_path = fullfile( projectroot, 'data','raw','exp', filesep );

% filenames of data to be processed
% full field measurements
list = {'101x101p_Chirp10_0-250kHz', ...                % 1                           
        '251x251p_0-1200kHz_Chirp_10us_x100', ...       % 2  
        '251x251p_400kHz_10HC_x10', ...                 % 3 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_200Hz', ...  % 4 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_5000Hz', ... % 5 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_500Hz', ...  % 6 
        '289x289p_0-350kHz_CHP160_x3_18Vpp_50Hz', ...   % 7 
        '289x289p_0-350kHz_CHP160_x3_2Vpp_200Hz', ...   % 8
        '289x289p_0-350kHz_CHP20_x3_18Vpp_50Hz', ...    % 9
        '289x289p_0-350kHz_CHP40_x3_18Vpp_50Hz', ...    % 10 
        '289x289p_0-350kHz_CHP80_x3_18Vpp_50Hz', ...    % 11 
        '289x289p_HANN500_x3_18Vpp_50Hz', ...           % 12 
        '289x289p_10HC400kHz_x3_7Vpp_50Hz', ...         % 13
        '289x289p_HANN100_x30_10Vpp_200Hz', ...         % 14 
        '289x289p_HANN100_x3_10Vpp_200Hz', ...          % 15 
        '289x289p_HANN25_x3_10Vpp_200Hz', ...           % 16 
        '289x289p_HANN50_x3_10Vpp_200Hz', ...           % 17 
        '289x289p_HANN50_x3_15Vpp_200Hz', ...           % 18 
        '289x289p_HANN50_x3_20Vpp_200Hz', ...           % 19 
        '289x289p_HANN50_x3_5Vpp_200Hz', ...            % 20 
        '493x493p_0-350kHz_CHP160_x3_18Vpp_50Hz', ...   % 21 
        '493x493p_HANN100_x10_10Vpp_200Hz', ...         % 22 
        '493z493p_0-350kHz_CHP20N_x3', ...              % 23
        '493z493p_0-350kHz_CHP40N_x3', ...              % 24 
        '493z493p_0-500kHz_CHP_40N_x3'};                % 25

fprintf('Making figures for: %s\n', foldername);
folder  = raw_data_path;
nFile   = length(test_case);
success = false(1, nFile);
for k = test_case
    filename = list{k};
    
    figfilename = [filename,'_KY_freq'];
    % check if already exist
    if(overwrite||(~overwrite && ~exist([output_path,figfilename,'.png'], 'file')))
        try 
            % load interim experimental data file
            disp('loading data');
            load([raw_data_path,filename]); % Data Length Width time
            nsY = 1024;
            nsX = 1024;
            nsT = 1024;
            %% PROCESS DATA
            fprintf('Processing:\n%s\n',filename);
            [KXKYF_,kx_vec,ky_vec,f_vec] = spatial_to_wavenumber_wavefield(Data,Length,Width,time);% only positive quarter
            %% KY-freq
            KYF = squeeze(abs(KXKYF_(:,1,:))); % @ kx=0
            wavenumber = ky_vec;
            if(~radians)
                wavenumber = wavenumber/(2*pi); % linear scale [1/m]
            end
            %% END OF PROCESSING
            %% START PLOTTING
            figure(k);
            set(gcf,'Color','w');
            imagesc(f_vec/1e3, wavenumber, KYF); set(gca,'YDir','normal');  
            KYFrow = reshape(KYF,(floor(nsY/2))*(floor(nsT/2)),1);
            KYFmax = max(KYFrow);
            x_hist = 0:KYFmax/10000:KYFmax;
            [n_elements, xout] = hist(KYFrow,x_hist);
            c_elements = cumsum(n_elements)/((floor(nsY/2))*(floor(nsT/2)));

            i = 1;
            while c_elements(i) < 0.999
                i = i+1;
            end 
            %b = max(caxis);  a = min(caxis); %b = 300;  
            axis([0 400 0 250]);
            %xlim([0 400]); ylim([0 250]);
            set(gca,'Fontsize',10,'linewidth',1);
            xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
            if(radians)
                ylabel({'$k_y$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
            else
                ylabel({'$k_y$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
            end
            colormap jet; 
            
            %caxis([0 max(caxis)/3]); % based on normalization to max
            caxis([0 xout(i)]); % based on histogram
            set(gca,'FontName','Times');
            fig = gcf;

            %title({[num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
            title({'3D FFT @ $k_x$ = 0'},'Fontsize',12,'interpreter','latex') 
            %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
            set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
            % remove unnecessary white space
            set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
            fig.PaperPositionMode   = 'auto';
            print([output_path,figfilename],'-dpng', '-r600'); 

            fprintf('Successfully processed:\n%s\n', filename);% successfully processed
             %% END PLOTTING
        catch
            fprintf('Failed: %s\n', filename);
        end
    else
        fprintf('Figure: %s already exist\n', figfilename);
    end
end
close all;



