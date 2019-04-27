% genetic algorithm - identification of elastic constants in fibre reinforced composite laminate

clear all; close all;

% load projectroot path
load project_paths projectroot src_path;
overwrite = false; % allow overwriting existing results if true
number_of_modes_considered = 4; % number of modes considered in calculation of objective function score
%% Load parameters which are used in experiment
% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );
% filename of parameter data
 filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_param'; 
 load([data_path,filename]); % wavenumber_max fmax beta number_of_wavenumber_points
% load experimental data file
exp_filename = {'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF'};
load([data_path,exp_filename{1}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points  
%% Input for SASE
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
layup = [0 90 0 90 90 0 90 0];
nlayers = length(layup);
h = [zeros(nlayers,1)+1]* 3e-3/nlayers; % thickness of layers;
% Stacking direction
stack_dir = 1;
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
[ObjV] = obj_ga_unidirectional(Phen,Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered);


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
       [ObjVSel] = obj_ga_unidirectional(bs2rv(SelCh,FieldD),Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered);
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
        plot(PBest(:,7));
        title('Volume fraction');
        figure(2);
        plot(Best);
        title(['Objectiv fun value, generation: ',num2str(gen)])
        drawnow;
        toc
   end 
   