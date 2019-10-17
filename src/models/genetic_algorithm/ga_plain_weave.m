% genetic algorithm - identification of elastic constants in fibre reinforced composite laminate

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true
number_of_modes_considered = 4; % number of modes considered in calculation of objective function score
%% Load parameters which are used in experiment
% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );
input_file = 2;
% filename of parameter data
 filename = {'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_param',...
                    'polar_interim_499x499p_chp200_x40_18Vpp_250Hz_KXKYF_param'}; 
 load([data_path,filename{input_file}]); % wavenumber_max fmax beta number_of_wavenumber_points
% load experimental data file
exp_filename = {'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF',... % 1 old plain weave
                            'polar_interim_499x499p_chp200_x40_18Vpp_250Hz_KXKYF'};         % 2 new plain weave
load([data_path,exp_filename{input_file}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points  
%% Input for SASE
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
%ht = 3/1000; % [m] laminate total thickness; old plain weave
ht = 3.9/1000; % [m] laminate total thickness; new plain weave
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
%layup = [0 0 0 0 0 0 0 0];
layup =[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
nlayers = length(layup);

h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers; new plain weave
% Stacking direction
stack_dir = 1;
%% input material properties for woven fabric
'weave type';                     fiberType = 'plainWeave';
'lamina thickness [mm]';     h_p = ht/nlayers; % 
'thickness of the fill [mm]';      h_f = h_p/2;
'thickness of the warp [mm]';      h_w = h_p/2;
'width of the fill [mm]';          a_f = 1.9;
'width of the warp [mm]';      a_w = 2;
'width of the fill gap [mm]';      g_f = 0.05;
'width of the warp gap [mm]';  g_w = 0.05;   
%% input for optimization
% lower and upper bounds of variables
%run(['inputs',filesep,'input1.m']);  % initial material properties
rhom0 = 1250; % kg/m^3
rhof0 = 1900; % kg/m^3
em0 = 3.43e9; % Pa
ef0 = 240e9; % Pa
nim0 = 0.35;
nif0 =  0.2; 
vol0 = 0.5;

variation = 0.2;
rhom_lb = (1-variation)*rhom0; % lower bound of matrix density
rhom_ub = (1+variation)*rhom0; % upper bound of matrix density
rhof_lb = (1-variation)*rhof0; % lower bound of fibres density
rhof_ub = (1+variation)*rhof0; % upper bound of fibres density
em_lb = (1-variation)*em0; % lower bound of Young's modulus of matrix
em_ub = (1+variation)*em0; % upper bound of Young's modulus of matrix
ef_lb = (1-variation)*ef0; % lower bound of Young's modulus of fibres
ef_ub = (1+variation)*ef0; % upper bound of Young's modulus of fibres
nim_lb = (1-variation)*nim0; % lower bound of Poisson's ratio of matrix
nim_ub = (1+variation)*nim0; % upper bound of Poisson's ratio of matrix
nif_lb = (1-variation)*nif0; % lower bound of Poisson's ratio of fibres
nif_ub = (1+variation)*nif0; % upper bound of Poisson's ratio of fibres
vol_lb = (1-variation)*vol0; % lower bound of volume fraction
vol_ub = (1+variation)*vol0; % upper bound of volume fraction
%% genetic algorithm parameters
NIND = 40;           % Number of individuals per subpopulations
MAXGEN = 40;        % maximum Number of generations
GGAP = 0.9;           % Generation gap, how many new individuals are created
NVAR = 7;           %number of variables in objective function
PRECI = 12;          % Precision of binary representation of variables

lb=[rhom_lb,rhof_lb,em_lb,ef_lb,nim_lb,nif_lb,vol_lb]; % lower bound for variables
ub=[rhom_ub,rhof_ub,em_ub,ef_ub,nim_ub,nif_ub,vol_ub]; % upper bound for variables
code=[1,1,1,1,1,1,1]; % Gray coding
scale=[0,0,0,0,0,0,0]; %arithmetic scale
lbin=  [1,1,1,1,1,1,1];%include lower bound of variable range
ubin= [1,1,1,1,1,1,1];%include upper bound of variable range
%%
%% tests loop
%%
for k_test = 1:1%100
    [k_test,100]
% Build field descriptor
%FieldD = [rep([PRECI],[1,NVAR]);lb;ub;code;scale;lbin;ubin];
FieldD = [rep([PRECI],[1,NVAR]);lb;ub; rep([1;0;1;1],[1,NVAR])];
% Initialise population (create random chromosomes)
Chrom = crtbp(NIND, NVAR*PRECI);
Phen = bs2rv(Chrom,FieldD); % convert binary to real
% Reset counters
   Best = NaN*ones(MAXGEN,1);	% best in current population
   PBest = NaN*ones(MAXGEN,7);	% best in current population
   gen = 0;			% generational counter

% Evaluate initial population
[ObjV] = obj_ga_plain_weave(Phen,Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,h_p,h_f,h_w,a_f,a_w,g_f,g_w,fiberType);

% Generational loop
   while gen < MAXGEN
        tic;
    % Assign fitness-value to entire population
       FitnV = ranking(ObjV);

    % Select individuals for breeding
       SelCh = select('sus', Chrom, FitnV, GGAP);

    % Recombine selected individuals (crossover)
       SelCh = recombin('xovsp',SelCh,0.7);

    % Perform mutation on offspring
       SelCh = mut(SelCh);

    % Evaluate offspring, call objective function
        %tic;
       [ObjVSel] = obj_ga_plain_weave(bs2rv(SelCh,FieldD),Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,h_p,h_f,h_w,a_f,a_w,g_f,g_w,fiberType);
        %toc
       % Reinsert offspring into current population
       [Chrom, ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);

    % Increment generational counter
       gen = gen+1;
       fprintf('Generation number: %d\n', gen);
    % Update display and record current best individual
        [Best(gen),I] = min(ObjV);
        fprintf('Best individual objective function value: %d\n', Best(gen));
        
        P=bs2rv(Chrom,FieldD);
        PBest(gen,:) = P(I,:);
        figure(1);
        plot(PBest(:,7),'o-');
        title('Volume fraction');
        figure(2);
        plot(Best,'o-');
        title(['Objective fun value, generation: ',num2str(gen)])
        drawnow;
        toc
   end 
   %save(['test_',num2str(k_test),'_3.9mm_16lay']);
%    save(['test_',num2str(k_test),'_3.9mm_8lay_Fabric_3']);
% save(['test_',num2str(k_test),'_4mm_8lay_Fabric_5']);
%save(['test_',num2str(k_test),'_3.9mm_8lay_plain_wave.mat']);

 %% Plot best case
radians = false;
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 8; 
[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);
load project_paths projectroot src_path;
%run([src_path,filesep,'models',filesep,'SASE',filesep,'inputs',filesep,'Fabric_6.m']);
rho_m = PBest(MAXGEN,1);
rho_f = PBest(MAXGEN,2);
e11_m = PBest(MAXGEN,3)/1e9;
e11_f = PBest(MAXGEN,4)/1e9;
ni12_m = PBest(MAXGEN,5);
ni12_f = PBest(MAXGEN,6);
vol_0 = PBest(MAXGEN,7);
e22_f = 0.1*e11_f;
ni23_f = PBest(MAXGEN,6);
format long;
[rho_m rho_f ]
[e11_m e11_f ]
[ni12_m ni12_f]
vol_0
ObjVal=min(ObjV);
%% Mechanical properties  
 [C11,C12,C13,C21,C22,C23,C31,C32,C33,C44,C55,C66,rho] = ...
        compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_0, ...
        e11_m, ni12_m, rho_m, e11_f, e22_f, ni12_f, ni23_f, rho_f,false);
 [C11,C12,C13,C22,C23,C33,C44,C55,C66]
 save(['out',filesep,'test_',num2str(k_test),'_3.9mm_',num2str(nlayers),'lay_plain_wave.mat'],'rho_m','rho_f','e11_m','e11_f','ni12_m','ni12_f','vol_0','C11','C12','C13','C22','C23','C33','C44','C55','C66','rho','ObjVal');
end

%% SASE
[wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);

for j=1:number_of_angles
    figure(j)
        set(gcf,'Color','w');
        %imagesc(fvec(2:end)/1e3, wavenumber(2:end,j)/2, squeeze(abs(Data_polar(j,2:end,2:end)))); 
        imagesc(fvec(2:end)/1e3, wavenumber(2:end,j), squeeze(abs(Data_polar(j,2:end,2:end)))); 
        set(gca,'YDir','normal'); 
        %axis([0 600 0 2000]);
        %axis([0 350 0 2000]);
        axis([0 350 0 min(wavenumber_max)]);
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
        if(radians)
            ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
        else
            ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
        end
        colormap jet; 
        %axis tight; 
        caxis([0 max(caxis)/3]); 

        hold on;
        fvec1=squeeze(FREQ(1,:,j)); % mode 1, angle j
        fvec2=squeeze(FREQ(2,:,j)); % mode 2, angle j
        fvec3=squeeze(FREQ(3,:,j)); % mode 3, angle j
        fvec4=squeeze(FREQ(4,:,j)); % mode 4, angle j
        kvec=squeeze(wavenumber(:,j)); % angle j
        plot(fvec1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec2(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec3(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec4(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        % get(gca,'FontName'); % default 'Helvetica'
        %set(gca,'FontName','Arial');
        %set(gca,'FontName','Helvetica');
        set(gca,'FontName','Times');
        fig = gcf;
        %title(['angle ', num2str(beta(j)),' deg'],'Fontsize',12);
        %title(['angle ', num2str(beta(j))],'Fontsize',12);
        %title(['angle ', num2str(beta(j)),' deg'],'Fontsize',12,'interpreter','latex');
        title({[num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
end  

%% Plot first generation/initial param
radians = false;
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 8; 
[number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
fvec = linspace(0,fmax,number_of_frequency_points);
load project_paths projectroot src_path;
%run([src_path,filesep,'models',filesep,'SASE',filesep,'inputs',filesep,'Fabric_5.m']);
% rho_m = PBest(1,1);
% rho_f = PBest(1,2);
% e11_m = PBest(1,3)/1e9;
% e11_f = PBest(1,4)/1e9;
% ni12_m = PBest(1,5);
% ni12_f = PBest(1,6);
% vol_0 = PBest(1,7);
% e22_f = 0.1*e11_f;
% ni23_f = PBest(1,6);
% rho_m = PBest(1,1);
rho_f = rhof0;
e11_m = em0/1e9;
e11_f = ef0/1e9;
ni12_m = nim0;
ni12_f = nif0;
vol_0 = vol0;
e22_f = 0.1*e11_f;
ni23_f = ni12_f;
%% Mechanical properties  
 [C11,C12,C13,C21,C22,C23,C31,C32,C33,C44,C55,C66,rho] = ...
        compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol_0, ...
        e11_m, ni12_m, rho_m, e11_f, e22_f, ni12_f, ni23_f, rho_f,false);
%% SASE
[wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);

for j=1:number_of_angles
    figure(number_of_angles+j)
        set(gcf,'Color','w');
        %imagesc(fvec(2:end)/1e3, wavenumber(2:end,j)/2, squeeze(abs(Data_polar(j,2:end,2:end)))); 
        imagesc(fvec(2:end)/1e3, wavenumber(2:end,j), squeeze(abs(Data_polar(j,2:end,2:end)))); 
        set(gca,'YDir','normal'); 
        %axis([0 600 0 2000]);
        %axis([0 350 0 2000]);
        axis([0 350 0 min(wavenumber_max)]);
        set(gca,'Fontsize',10,'linewidth',1);
        xlabel({'$f$ [kHz]'},'Fontsize',12,'interpreter','latex');
        if(radians)
            ylabel({'$k$ [rad/m]'},'Fontsize',12,'interpreter','latex'); % radian scale [rad/m]
        else
            ylabel({'$k$ [1/m]'},'Fontsize',12,'interpreter','latex'); % linear scale [1/m]
        end
        colormap jet; 
        %axis tight; 
        caxis([0 max(caxis)/3]); 

        hold on;
        fvec1=squeeze(FREQ(1,:,j)); % mode 1, angle j
        fvec2=squeeze(FREQ(2,:,j)); % mode 2, angle j
        fvec3=squeeze(FREQ(3,:,j)); % mode 3, angle j
        fvec4=squeeze(FREQ(4,:,j)); % mode 4, angle j
        kvec=squeeze(wavenumber(:,j)); % angle j
        plot(fvec1(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec2(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec3(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        plot(fvec4(2:end)/1e3,kvec(2:end),'y','linewidth',1);
        % get(gca,'FontName'); % default 'Helvetica'
        %set(gca,'FontName','Arial');
        %set(gca,'FontName','Helvetica');
        set(gca,'FontName','Times');
        fig = gcf;
        %title(['angle ', num2str(beta(j)),' deg'],'Fontsize',12);
        %title(['angle ', num2str(beta(j))],'Fontsize',12);
        %title(['angle ', num2str(beta(j)),' deg'],'Fontsize',12,'interpreter','latex');
        title({[num2str(beta(j)),'$^{\circ}$']},'Fontsize',12,'interpreter','latex');
        %set(gca, 'Position',[0 0 1.2 1.2]); % figure without axis and white border
        set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % size 12cm by 8cm (1-column text)
        % remove unnecessary white space
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
        fig.PaperPositionMode   = 'auto';
end  