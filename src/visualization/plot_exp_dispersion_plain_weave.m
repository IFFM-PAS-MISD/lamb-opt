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
foldername = 'composite_structures_GA';
radians = false;
% create output path
output_path = prepare_exp_figure_paths(foldername);
% length and width of full wavefield area
Length = 0.726;    
Width = 0.726;
test_case=[1]; % select file numbers for processing

% create path to the experimental raw data folder
raw_data_path = fullfile( projectroot, 'data','raw','exp', filesep );

% filenames of data to be processed
% full field measurements
list = {'499x499p_chp200_x40_18Vpp_250Hz', ...          % 1                           
        '499x499p_chp200_x8_18Vpp_100Hz' };                % 2

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



