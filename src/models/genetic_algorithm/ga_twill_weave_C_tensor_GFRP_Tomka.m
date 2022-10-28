% genetic algorithm - identification of elastic constants in fibre reinforced composite laminate
% GFRP Tomka, twill weave

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;
%% Load parameters which are used in experiment
% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );

%% Prepare output directories
% allow overwriting existing results if true
overwrite=true;
% retrieve model name based on running file and folder
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
% prepare model output path
model_output_path = prepare_model_paths('raw','num',modelfolder,modelname);
figure_output_path = prepare_figure_paths(modelname);
% test_case 1 and 2: 4 modes considered
% other cases: 6 modes considered
number_of_modes_considered = 6; % number of modes considered in calculation of objective function score

%% Load parameters which are used in experiment
% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );
input_file = 1;
% filename of parameter data
 filename = {'polar_interim_chirp_interpLaminat_7_GFRP_Tomka_KXKYF_param'}; 
 load([data_path,filename{input_file}]); % wavenumber_max fmax beta number_of_wavenumber_points
% load experimental data file
exp_filename = {'polar_interim_chirp_interpLaminat_7_GFRP_Tomka_KXKYF'};     
load([data_path,exp_filename{input_file}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points  
%% Input for SASE
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
ht = 2/1000; % [m] laminate total thickness; 
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
layup = [0 90 0 90 0 90 90 0 90 0 90 0];
nlayers = length(layup);

h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers; new plain weave
% Stacking direction
stack_dir = 1;
%%
% known parameters
% m=8.55; % total mass of the specimen [kg]
% V=0.5*0.5*ht; % specimen volume [m^3]
rho = 1700; % [kg/m^3]
%% input for optimization
% lower and upper bounds of variables
Q11_0 = 26e9;
Q12_0 = 3.7e9;
Q13_0 = 3.7e9;
Q22_0 = 26e9;
Q23_0 = 3e9;
Q33_0 = 9.4e9;
Q44_0 = 3e9;
Q55_0 = 3e9;
Q66_0 = 3e9;

variation = 0.7;
Q11_lb = (1-variation)*Q11_0; 
Q11_ub = (1+variation)*Q11_0; 
Q12_lb = (1-variation)*Q12_0; 
Q12_ub = (1+variation)*Q12_0; 
Q13_lb = (1-variation)*Q13_0; 
Q13_ub = (1+variation)*Q13_0; 
Q22_lb = (1-variation)*Q22_0; 
Q22_ub = (1+variation)*Q22_0; 
Q23_lb = (1-variation)*Q23_0; 
Q23_ub = (1+variation)*Q23_0; 
Q33_lb = (1-variation)*Q33_0; 
Q33_ub = (1+variation)*Q33_0; 
Q44_lb = (1-variation)*Q44_0; 
Q44_ub = (1+variation)*Q44_0; 
Q55_lb = (1-variation)*Q55_0; 
Q55_ub = (1+variation)*Q55_0; 
Q66_lb = (1-variation)*Q66_0; 
Q66_ub = (1+variation)*Q66_0; 
%% genetic algorithm parameters
% fittnes function scaling factors
a=100;
b=310;
NIND = 100;           % Number of individuals per subpopulations
MAXGEN = 500;        % maximum Number of generations
GGAP = 0.9;           % Generation gap, how many new individuals are created
NVAR = 9;           %number of variables in objective function
PRECI = 12;          % Precision of binary representation of variables
XOV = 0.7;          % crossover rate
MUTR = 0.05;         % Mutation rate
lb=[Q11_lb,Q12_lb,Q13_lb,Q22_lb,Q23_lb,Q33_lb,Q44_lb,Q55_lb,Q66_lb]; % lower bound for variables
ub=[Q11_ub,Q12_ub,Q13_ub,Q22_ub,Q23_ub,Q33_ub,Q44_ub,Q55_ub,Q66_ub]; % upper bound for variables
code=[1,1,1,1,1,1,1,1,1]; % Gray coding
scale=[0,0,0,0,0,0,0,0,0]; %arithmetic scale
lbin=  [1,1,1,1,1,1,1,1,1];%include lower bound of variable range
ubin= [1,1,1,1,1,1,1,1,1];%include upper bound of variable range
%%
%% tests loop
%%
for test_case = [3:10]
    
    output_name = [model_output_path,filesep,num2str(test_case),'output'];
     if(overwrite||(~overwrite && ~exist([output_name,'.mat'], 'file')))
        fprintf([modelname,' test case: %d\n'], test_case);
        % Build field descriptor
        %FieldD = [rep([PRECI],[1,NVAR]);lb;ub;code;scale;lbin;ubin];
        FieldD = [rep([PRECI],[1,NVAR]);lb;ub; rep([1;0;1;1],[1,NVAR])];
        % Initialise population (create random chromosomes)
        Chrom = crtbp(NIND, NVAR*PRECI);
        Phen = bs2rv(Chrom,FieldD); % convert binary to real
        % Reset counters
       Best = NaN*ones(MAXGEN,1);	% best in current population
       Mean = NaN*ones(MAXGEN,1);	% mean in current population
       PBest = NaN*ones(MAXGEN,NVAR);	% best in current population
       ObjV_limit = 8.5; % stopping criteria
       gen = 0;			% generational counter

        % Evaluate initial population
        [ObjV] = obj_ga_C_tensor_known_mass(Phen,Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,rho,a,b);
        % Generational loop
       while gen < MAXGEN
            tic;
        % Assign fitness-value to entire population
           FitnV = ranking(ObjV);

        % Select individuals for breeding
           SelCh = select('sus', Chrom, FitnV, GGAP); % stochastic universal sampling
           %SelCh = select('rws', Chrom, FitnV, GGAP); % stochastic sampling with replacement
        % Recombine selected individuals (crossover)
           SelCh = recombin('xovsp',SelCh,0.7);

        % Perform mutation on offspring
           %SelCh = mut(SelCh); % Probability of mutation Pm=0.7/Lind where Lind is the length of chromosome structure
           SelCh = mut(SelCh,MUTR); % Probability of mutation MUTR
        % Evaluate offspring, call objective function
            %tic;
           [ObjVSel] = obj_ga_C_tensor_known_mass(bs2rv(SelCh,FieldD),Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,rho,a,b);
           %toc
           % Reinsert offspring into current population
           [Chrom, ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);

        % Increment generational counter
           gen = gen+1;
           fprintf('Generation number: %d\n', gen);
        % Update display and record current best individual
            [Best(gen),I] = min(ObjV);
            Mean(gen) = mean(ObjV);
            fprintf('Best individual objective function value: %d\n', Best(gen));

            P=bs2rv(Chrom,FieldD);
            PBest(gen,:) = P(I,:);
            close all;
            figure(1);
            subplot(3,3,1);
            plot(PBest(:,1),'o-');
            ylim([Q11_lb Q11_ub]);
            title('C11');
            subplot(3,3,2);
            plot(PBest(:,2),'o-');
            ylim([Q12_lb Q12_ub]);
            title('C12');
            subplot(3,3,3);
            plot(PBest(:,3),'o-');
            ylim([Q13_lb Q13_ub]);
            title('C13');
            subplot(3,3,4);
            plot(PBest(:,4),'o-');
            title('C22');
            subplot(3,3,5);
            plot(PBest(:,5),'o-');
            ylim([Q23_lb Q23_ub]);
            title('C23');
            subplot(3,3,6);
            plot(PBest(:,6),'o-');
            ylim([Q33_lb Q33_ub]);
            title('C33');
            subplot(3,3,7);
            plot(PBest(:,7),'o-');
            ylim([Q44_lb Q44_ub]);
            title('C44');
            subplot(3,3,8);
            plot(PBest(:,8),'o-');
            ylim([Q55_lb Q55_ub]);
            title('C55');
            subplot(3,3,9);
            plot(PBest(:,9),'o-');
            ylim([Q66_lb Q66_ub]);
            title('C66');
            figure(2);
            plot(Best,'bo-','MarkerSize',2);hold on;
            plot(Mean,'rd-','MarkerSize',2);
            legend('Best','Mean');
            title(['Objective fun value, generation: ',num2str(gen)])
            drawnow;
            toc
            if(Best(gen) < ObjV_limit && gen >= 50) 
                break; 
            end
       end 

         
        
        C11 = PBest(gen,1);
        C12 = PBest(gen,2);
        C13 = PBest(gen,3);
        C22 = PBest(gen,4);
        C23 = PBest(gen,5);
        C33 = PBest(gen,6);
        C44 = PBest(gen,7);
        C55 = PBest(gen,8);
        C66 = PBest(gen,9);
        ObjVal=min(ObjV);

        [C11,C12,C13,C22,C23,C33,C44,C55,C66]

         save(output_name,'C11','C12','C13','C22','C23','C33','C44','C55','C66','rho','ObjVal');
         %% Plot best case
         radians = false;
         % size 12cm by 8cm (1-column text)
         %fig_width = 12; fig_height = 8; 
          % size 7cm by 5cm (2-column text)
         fig_width = 7; fig_height = 5; 
         [number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
         fvec = linspace(0,fmax,number_of_frequency_points);
         [wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);
         if(~radians)
            wavenumber = wavenumber/(2*pi); % linear scale [1/m]
            wavenumber_max = wavenumber_max/(2*pi); % linear scale [1/m]
         end
         figure(3); % angle 1 0deg
         set(gcf,'Color','w');
         imagesc(fvec(2:end)/1e3, wavenumber(2:end,1), squeeze(abs(Data_polar(1,2:end,2:end)))); 
         set(gca,'YDir','normal'); 
         axis([0 500 0 400]);
         set(gca,'Fontsize',10,'linewidth',1);
         xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
         if(radians)
             ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
         else
             ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
         end
         %colormap jet; 
         colormap turbo; 
         %axis tight; 
         caxis([0 max(caxis)/3]); 

         hold on;
         fvec1=squeeze(FREQ(1,:,1)); % mode 1, angle j
         fvec2=squeeze(FREQ(2,:,1)); % mode 2, angle j
         fvec3=squeeze(FREQ(3,:,1)); % mode 3, angle j
         fvec4=squeeze(FREQ(4,:,1)); % mode 4, angle j
         fvec5=squeeze(FREQ(5,:,1)); % mode 5, angle j
         fvec6=squeeze(FREQ(6,:,1)); % mode 5, angle j
         kvec=squeeze(wavenumber(:,1)); % angle j
         LW=0.5; % small figures
         %LW=1; % large figures
         plot(fvec1(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % A0
         plot(fvec2(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % SH0
         plot(fvec3(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % S0
         plot(fvec4(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(fvec5(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(fvec6(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         %set(gca,'FontName','Helvetica');% default 'Helvetica'
         set(gca,'FontName','Times');
         fig = gcf;
         title({['$Gen=$ ',num2str(gen),', ',num2str(beta(1)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
         %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
         set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
         % remove unnecessary white space
         set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
         fig.PaperPositionMode   = 'auto';
         figfilename = ['Run_',num2str(test_case),'_Gen_',num2str(gen),'_angle_',num2str(beta(1))];
         print([figure_output_path,figfilename],'-dpng', '-r600'); 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         figure(4);
         set(gcf,'Color','w');
        t=tiledlayout(1,2);
        %t.TileSpacing = 'tight';
        t.TileSpacing = 'none';
        t.Padding = 'tight';
        
        % Left plot
        ax1 = nexttile;
        imagesc(ax1,fvec(2:end)/1e3, wavenumber(2:end,1), squeeze(abs(Data_polar(1,2:end,2:end)))); 
         set(gca,'YDir','normal'); 
         axis([0 500 0 400]);
         set(gca,'Fontsize',10,'linewidth',1);
         xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
         if(radians)
             ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
         else
             ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
         end
         %colormap jet; 
         colormap turbo;
         %axis tight; 
         caxis([0 max(caxis)/3]); 

         hold on;
         fvec1=squeeze(FREQ(1,:,1)); % mode 1, angle j
         fvec2=squeeze(FREQ(2,:,1)); % mode 2, angle j
         fvec3=squeeze(FREQ(3,:,1)); % mode 3, angle j
         fvec4=squeeze(FREQ(4,:,1)); % mode 4, angle j
         fvec5=squeeze(FREQ(5,:,1)); % mode 5, angle j
         fvec6=squeeze(FREQ(6,:,1)); % mode 5, angle j
         kvec=squeeze(wavenumber(:,1)); % angle j
         LW=0.5; % small figures
         %LW=1; % large figures
         plot(ax1,fvec1(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % A0
         plot(ax1,fvec2(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % SH0
         plot(ax1,fvec3(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % S0
         plot(ax1,fvec4(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax1,fvec5(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax1,fvec6(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         %set(gca,'FontName','Helvetica');% default 'Helvetica'
         set(ax1,'FontName','Times');
         title({[num2str(beta(1)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
         text(ax1,0,440,['Gen= ',num2str(gen)],'HorizontalAlignment','left','Fontsize',10,'FontName','Times','fontweight','bold');
            
         % Right plot
         ax2 = nexttile;
         imagesc(ax2,fvec(2:end)/1e3, wavenumber(2:end,4), squeeze(abs(Data_polar(4,2:end,2:end)))); 
         set(ax2,'YDir','normal'); 
         axis([0 500 0 400]);
         set(ax2,'Fontsize',10,'linewidth',1);
         xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
%          if(radians)
%              ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
%          else
%              ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
%          end
         %colormap jet; 
         colormap turbo;
         %axis tight; 
         caxis([0 max(caxis)/3]); 

         hold on;
         fvec1=squeeze(FREQ(1,:,4)); % mode 1, angle 4
         fvec2=squeeze(FREQ(2,:,4)); % mode 2, angle 4
         fvec3=squeeze(FREQ(3,:,4)); % mode 3, angle 4
         fvec4=squeeze(FREQ(4,:,4)); % mode 4, angle 4
         fvec5=squeeze(FREQ(5,:,4)); % mode 5, angle 4
         fvec6=squeeze(FREQ(6,:,4)); % mode 5, angle 4
         kvec=squeeze(wavenumber(:,4)); % angle 4
         LW=0.5; % small figures
         %LW=1; % large figures
         plot(ax2,fvec1(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % A0
         plot(ax2,fvec2(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % SH0
         plot(ax2,fvec3(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % S0
         plot(ax2,fvec4(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax2,fvec5(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax2,fvec6(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         %set(gca,'FontName','Helvetica');% default 'Helvetica'
         set(ax2,'FontName','Times');
         title({[num2str(beta(4)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
         fig = gcf;
         
         set(fig, 'Units','centimeters', 'Position',[10 10 2*fig_width fig_height]); % 
         % remove unnecessary white space
         set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
         fig.PaperPositionMode   = 'auto';
         figfilename = ['Run_',num2str(test_case),'_Gen_',num2str(gen),'_2angles'];
         print([figure_output_path,figfilename],'-dpng', '-r600'); 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         figure(5);
         set(gcf,'Color','w');
        t=tiledlayout(2,2);
        %t.TileSpacing = 'tight';
        t.TileSpacing = 'none';
        t.Padding = 'tight';
        
        % plot (1,1) 0deg
        ax1 = nexttile;
        imagesc(ax1,fvec(2:end)/1e3, wavenumber(2:end,1), squeeze(abs(Data_polar(1,2:end,2:end)))); 
         set(gca,'YDir','normal'); 
         axis([0 500 0 400]);
         set(gca,'Fontsize',10,'linewidth',1);
         %xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
         if(radians)
             ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
         else
             ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
         end
         %colormap jet; 
         colormap turbo;
         %axis tight; 
         caxis([0 max(caxis)/3]); 

         hold on;
         fvec1=squeeze(FREQ(1,:,1)); % mode 1, angle j
         fvec2=squeeze(FREQ(2,:,1)); % mode 2, angle j
         fvec3=squeeze(FREQ(3,:,1)); % mode 3, angle j
         fvec4=squeeze(FREQ(4,:,1)); % mode 4, angle j
         fvec5=squeeze(FREQ(5,:,1)); % mode 5, angle j
         fvec6=squeeze(FREQ(6,:,1)); % mode 5, angle j
         kvec=squeeze(wavenumber(:,1)); % angle j
         LW=0.5; % small figures
         %LW=1; % large figures
         plot(ax1,fvec1(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % A0
         plot(ax1,fvec2(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % SH0
         plot(ax1,fvec3(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % S0
         plot(ax1,fvec4(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax1,fvec5(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax1,fvec6(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         %set(gca,'FontName','Helvetica');% default 'Helvetica'
         set(ax1,'FontName','Times');
         title({[num2str(beta(1)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
         text(ax1,0,440,['Gen= ',num2str(gen)],'HorizontalAlignment','left','Fontsize',10,'FontName','Times','fontweight','bold');
            
         % plot (1,2) 30 deg
         ax2 = nexttile;
         imagesc(ax2,fvec(2:end)/1e3, wavenumber(2:end,3), squeeze(abs(Data_polar(3,2:end,2:end)))); 
         set(ax2,'YDir','normal'); 
         axis([0 500 0 400]);
         set(ax2,'Fontsize',10,'linewidth',1);
         %xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
%          if(radians)
%              ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
%          else
%              ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
%          end
         %colormap jet; 
         colormap turbo;
         %axis tight; 
         caxis([0 max(caxis)/3]); 

         hold on;
         fvec1=squeeze(FREQ(1,:,3)); % mode 1, angle 3
         fvec2=squeeze(FREQ(2,:,3)); % mode 2, angle 3
         fvec3=squeeze(FREQ(3,:,3)); % mode 3, angle 3
         fvec4=squeeze(FREQ(4,:,3)); % mode 4, angle 3
         fvec5=squeeze(FREQ(5,:,3)); % mode 5, angle 3
         fvec6=squeeze(FREQ(6,:,3)); % mode 5, angle 3
         kvec=squeeze(wavenumber(:,3)); % angle 3
         LW=0.5; % small figures
         %LW=1; % large figures
         plot(ax2,fvec1(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % A0
         plot(ax2,fvec2(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % SH0
         plot(ax2,fvec3(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % S0
         plot(ax2,fvec4(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax2,fvec5(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax2,fvec6(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         %set(gca,'FontName','Helvetica');% default 'Helvetica'
         set(ax2,'FontName','Times');
         title({[num2str(beta(3)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
         
         % plot (2,1) 60deg
        ax3 = nexttile;
        imagesc(ax3,fvec(2:end)/1e3, wavenumber(2:end,5), squeeze(abs(Data_polar(5,2:end,2:end)))); 
         set(gca,'YDir','normal'); 
         axis([0 500 0 400]);
         set(gca,'Fontsize',10,'linewidth',1);
         xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
         if(radians)
             ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
         else
             ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
         end
         %colormap jet; 
         colormap turbo;
         %axis tight; 
         caxis([0 max(caxis)/3]); 

         hold on;
         fvec1=squeeze(FREQ(1,:,5)); % mode 1, angle 5
         fvec2=squeeze(FREQ(2,:,5)); % mode 2, angle 5
         fvec3=squeeze(FREQ(3,:,5)); % mode 3, angle 5
         fvec4=squeeze(FREQ(4,:,5)); % mode 4, angle 5
         fvec5=squeeze(FREQ(5,:,5)); % mode 5, angle 5
         fvec6=squeeze(FREQ(6,:,5)); % mode 5, angle 5
         kvec=squeeze(wavenumber(:,5)); % angle j
         LW=0.5; % small figures
         %LW=1; % large figures
         plot(ax3,fvec1(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % A0
         plot(ax3,fvec2(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % SH0
         plot(ax3,fvec3(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % S0
         plot(ax3,fvec4(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax3,fvec5(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax3,fvec6(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         %set(gca,'FontName','Helvetica');% default 'Helvetica'
         set(ax3,'FontName','Times');
         title({[num2str(beta(5)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
            
         % plot (2,2) 90 deg
         ax4 = nexttile;
         imagesc(ax4,fvec(2:end)/1e3, wavenumber(2:end,7), squeeze(abs(Data_polar(7,2:end,2:end)))); 
         set(ax4,'YDir','normal'); 
         axis([0 500 0 400]);
         set(ax4,'Fontsize',10,'linewidth',1);
         xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
%          if(radians)
%              ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
%          else
%              ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
%          end
         %colormap jet; 
         colormap turbo;
         %axis tight; 
         caxis([0 max(caxis)/3]); 

         hold on;
         fvec1=squeeze(FREQ(1,:,7)); % mode 1, angle 7
         fvec2=squeeze(FREQ(2,:,7)); % mode 2, angle 7
         fvec3=squeeze(FREQ(3,:,7)); % mode 3, angle 7
         fvec4=squeeze(FREQ(4,:,7)); % mode 4, angle 7
         fvec5=squeeze(FREQ(5,:,7)); % mode 5, angle 7
         fvec6=squeeze(FREQ(6,:,7)); % mode 5, angle 7
         kvec=squeeze(wavenumber(:,7)); % angle 7
         LW=0.5; % small figures
         %LW=1; % large figures
         plot(ax4,fvec1(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % A0
         plot(ax4,fvec2(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % SH0
         plot(ax4,fvec3(2:end)/1e3,kvec(2:end),'w','linewidth',LW); % S0
         plot(ax4,fvec4(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax4,fvec5(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         plot(ax4,fvec6(2:end)/1e3,kvec(2:end),'w','linewidth',LW);
         %set(gca,'FontName','Helvetica');% default 'Helvetica'
         set(ax4,'FontName','Times');
         title({[num2str(beta(7)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
         
         fig = gcf;
         set(fig, 'Units','centimeters', 'Position',[10 10 2*fig_width 2*fig_height]); % 
         % remove unnecessary white space
         set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
         fig.PaperPositionMode   = 'auto';
         figfilename = ['Run_',num2str(test_case),'_Gen_',num2str(gen),'_4angles'];
         print([figure_output_path,figfilename],'-dpng', '-r600'); 
     else
         fprintf([modelname,' test case: %d already exist\n'], test_case);
     end
     %plot GA convergence  
    figure(6);
    plot(Best,'bo-','MarkerFaceColor','b','MarkerSize',2);hold on;
    plot(Mean,'rd-','MarkerFaceColor','r','MarkerSize',2);
    legend('Best','Mean');
    
    set(gcf,'Color','w');
    set(gca,'Fontsize',10,'linewidth',1);
    title('');
    xlabel('Gen','Fontsize',12);
    ylabel('F','Fontsize',12);
    %axis([0 50 -10 100]);
    set(gca,'FontName','Times');
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    figfilename = ['Run_',num2str(test_case),'_GA_convergence'];
    print([paper_path,figfilename],'-dpng', '-r600');  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(7);
    set(gcf,'Color','w');
    t=tiledlayout(3,3);
    %t.TileSpacing = 'tight';
    t.TileSpacing = 'none';
    t.Padding = 'tight';
    
    ax1 = nexttile;
    plot(ax1,PBest(:,1)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q11_lb/1e9 Q11_ub/1e9]);
    title({'$C_{11}$'},'Fontsize',12,'interpreter','latex');
    ax2 = nexttile;
    plot(ax2,PBest(:,2)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q12_lb/1e9 Q12_ub/1e9]);
    title({'$C_{12}$'},'Fontsize',12,'interpreter','latex');
    ax3 = nexttile;
    plot(ax3,PBest(:,3)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q13_lb/1e9 Q13_ub/1e9]);
    title({'$C_{13}$'},'Fontsize',12,'interpreter','latex');
    ax4 = nexttile;
    plot(ax4,PBest(:,4)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q22_lb/1e9 Q22_ub/1e9]);
    title({'$C_{22}$'},'Fontsize',12,'interpreter','latex');
    ax5 = nexttile;
    plot(ax5,PBest(:,5)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q23_lb/1e9 Q23_ub/1e9]);
    title({'$C_{23}$'},'Fontsize',12,'interpreter','latex');
    ax6 = nexttile;
    plot(ax6,PBest(:,6)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q33_lb/1e9 Q33_ub/1e9]);
    title({'$C_{33}$'},'Fontsize',12,'interpreter','latex');
    ax7 = nexttile;
    plot(ax7,PBest(:,7)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q44_lb/1e9 Q44_ub/1e9]);
    title({'$C_{44}$'},'Fontsize',12,'interpreter','latex');
    ax8 = nexttile;
    plot(ax8,PBest(:,8)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q55_lb/1e9 Q55_ub/1e9]);
    title({'$C_{55}$'},'Fontsize',12,'interpreter','latex');
    ax9 = nexttile;
    plot(ax9,PBest(:,9)/1e9,'o-','MarkerFaceColor','k','MarkerSize',2);
    ylim([Q66_lb/1e9 Q66_ub/1e9]);
    title({'$C_{66}$'},'Fontsize',12,'interpreter','latex');
    set(gca,'FontName','Times');
    fig = gcf;
    set(fig, 'Units','centimeters', 'Position',[10 10 2*fig_width 2*fig_height]); % 
    % remove unnecessary white space
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    fig.PaperPositionMode   = 'auto';
    figfilename = ['Run_',num2str(test_case),'_C_convergence'];
    print([paper_path,figfilename],'-dpng', '-r600');  
            
end
