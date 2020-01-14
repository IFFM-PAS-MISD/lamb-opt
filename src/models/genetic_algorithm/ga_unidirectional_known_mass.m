% genetic algorithm - identification of elastic constants in fibre reinforced composite laminate

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
number_of_modes_considered = 4; % number of modes considered in calculation of objective function score
%% Load parameters which are used in experiment
% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );
input_file = 1;
% filename of parameter data
 filename = {'polar_interim_499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni_KXKYF_param',...
                    'polar_interim_499x499p_chp200_x40_6Vpp_250Hz_uni_KXKYF_param'}; 
 load([data_path,filename{input_file}]); % wavenumber_max fmax beta number_of_wavenumber_points
% load experimental data file
exp_filename = {'polar_interim_499x499p_chp200_x30_6Vpp_250Hz_100mmsv_small_uni_KXKYF',... % 1 small area unidirectional
                            'polar_interim_499x499p_chp200_x40_6Vpp_250Hz_uni_KXKYF'};         % 2 large area unidirectional
load([data_path,exp_filename{input_file}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points  
%% Input for SASE
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
ht = 2.85/1000; % [m] laminate total thickness; unidirectional
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
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
%% input for optimization
% lower and upper bounds of variables
%run(['inputs',filesep,'input1.m']);  % initial material properties
rhof0 = 1900; % kg/m^3
em0 = 3.43e9; % Pa
ef0 = 240e9; % Pa
nim0 = 0.35;
nif0 =  0.2; 
%vol0 = 0.5;
vol0 = 0.45;

variation = 0.25;
rhom_lb = 1150; % lower bound of cured resing density (from manufacturer)
rhom_ub = 1250; % upper bound of cured resing density (from manufacturer)
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
NIND = 100;           % Number of individuals per subpopulations
MAXGEN = 70;        % maximum Number of generations
GGAP = 0.9;           % Generation gap, how many new individuals are created
NVAR = 6;           %number of variables in objective function
PRECI = 12;          % Precision of binary representation of variables

lb=[rhom_lb,em_lb,ef_lb,nim_lb,nif_lb,vol_lb]; % lower bound for variables
ub=[rhom_ub,em_ub,ef_ub,nim_ub,nif_ub,vol_ub]; % upper bound for variables
code=[1,1,1,1,1,1]; % Gray coding
scale=[0,0,0,0,0,0]; %arithmetic scale
lbin=  [1,1,1,1,1,1];%include lower bound of variable range
ubin= [1,1,1,1,1,1];%include upper bound of variable range
%%
%% tests loop
%%
for test_case = 4%[1:50]
    
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
       ObjV_limit = 8.0; % stopping criteria
       gen = 0;			% generational counter

        % Evaluate initial population
        [ObjV] = obj_ga_unidirectional_known_mass(Phen,Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,rho);

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
           SelCh = mut(SelCh);

        % Evaluate offspring, call objective function
            %tic;
           [ObjVSel] = obj_ga_unidirectional_known_mass(bs2rv(SelCh,FieldD),Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,rho);

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
            figure(1);
            plot(PBest(:,6),'o-');
            title('Volume fraction');
            figure(2);
            plot(Best,'bo-');hold on;
            plot(Mean,'rd-');
            legend('Best','Mean');
            title(['Objective fun value, generation: ',num2str(gen)])
            drawnow;
            toc
            if(Best(gen) < ObjV_limit && gen >= 50) 
                break; 
            end
       end 

         %% Plot best case
        radians = false;
        % size 12cm by 8cm (1-column text)
        fig_width = 12; fig_height = 8; 
        [number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
        fvec = linspace(0,fmax,number_of_frequency_points);
        load project_paths projectroot src_path;
       
        rhom =  PBest(gen,1);
        em =  PBest(gen,2);
        ef = PBest(gen,3);
        nim = PBest(gen,4);
        nif = PBest(gen,5);
        vol = PBest(gen,6);
        rhof = (rho - rhom*(1-vol))/vol;     
        format long;
        ObjVal=min(ObjV);
        %% Mechanical properties  
         % homogenization of unidirectional fibre reinforce composite
        [rho,e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23] = homogenization(rhom,rhof,em,ef,nim,nif,vol);
        % Elastic constants of composite lamina in terms of principal material directions
        [C11,C12,C13,C22,C23,C33,C44,C55,C66] = lamina_3D(e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23);   
         save(output_name,'rhom','rhof','em','ef','nim','nif','vol','C11','C12','C13','C22','C23','C33','C44','C55','C66','rho','ObjVal');
     else
         fprintf([modelname,' test case: %d already exist\n'], test_case);
     end
end

 