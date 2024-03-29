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
number_of_modes_considered = 2; % number of modes considered in calculation of objective function score
%% Load parameters which are used in experiment
% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );
input_file = 1;
% filename of parameter data
 filename = {'polar_interim_483x483p_CHIRP_20-500kHz_125us_6Vpp_x3_stringer_intact_KXKYF_param'}; 
 load([data_path,filename{input_file}]); % wavenumber_max fmax beta number_of_wavenumber_points
% load experimental data file
exp_filename = {'polar_interim_483x483p_CHIRP_20-500kHz_125us_6Vpp_x3_stringer_intact_KXKYF'};
load([data_path,exp_filename{input_file}]); % Data_polar wavenumber_max fmax beta number_of_wavenumber_points  
%return;
%% Input for SASE
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
ht = 2.2/1000; % [m] laminate total thickness;
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
layup = [45,0,-45,90,-45,0,45,90,90,45,0,-45,90,-45,0,45];

nlayers = length(layup);

h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers; new plain weave
% Stacking direction
stack_dir = 1;
%%
% known parameters
m=6.46; % total mass of the specimen [kg]
V=1.2*1.2*ht; % specimen volume [m^3]
rho = 1571; % [kg/m^3]
%% input for optimization
% lower and upper bounds of variables
Q11_0 = 130e9; % 
Q12_0 = 6.0e9;
Q13_0 = 6.0e9;
Q22_0 = 11e9;
Q23_0 = 5e9;
Q33_0 = 11e9;
Q44_0 = 3e9;
Q55_0 = 4.2e9;
Q66_0 = 4.2e9;

variation = 0.5; % 20% is too small - change to 50%
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
NIND = 100;           % Number of individuals per subpopulations
MAXGEN = 120;        % maximum Number of generations
GGAP = 0.8;           % Generation gap, how many new individuals are created
NVAR = 9;           %number of variables in objective function
PRECI = 12;          % Precision of binary representation of variables

lb=[Q11_lb,Q12_lb,Q13_lb,Q22_lb,Q23_lb,Q33_lb,Q44_lb,Q55_lb,Q66_lb]; % lower bound for variables
ub=[Q11_ub,Q12_ub,Q13_ub,Q22_ub,Q23_ub,Q33_ub,Q44_ub,Q55_ub,Q66_ub]; % upper bound for variables
code=[1,1,1,1,1,1,1,1,1]; % Gray coding
scale=[0,0,0,0,0,0,0,0,0]; %arithmetic scale
lbin=  [1,1,1,1,1,1,1,1,1];%include lower bound of variable range
ubin= [1,1,1,1,1,1,1,1,1];%include upper bound of variable range
%%
%% tests loop
%%
for test_case = [1:5]
    
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
        [ObjV] = obj_ga_C_tensor_known_mass_unidirectional(Phen,Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,rho);
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
           SelCh = mut(SelCh); % Probability of mutation Pm=0.7/Lind where Lind is the length of chromosome structure

        % Evaluate offspring, call objective function
            %tic;
           [ObjVSel] = obj_ga_C_tensor_known_mass_unidirectional(bs2rv(SelCh,FieldD),Data_polar,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer,fmax,number_of_modes_considered,rho);
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
            subplot(3,3,1);
            plot(PBest(:,1),'o-');
            title('C11');
            subplot(3,3,2);
            plot(PBest(:,2),'o-');
            title('C12');
            subplot(3,3,3);
            plot(PBest(:,3),'o-');
            title('C13');
            subplot(3,3,4);
            plot(PBest(:,4),'o-');
            title('C22');
            subplot(3,3,5);
            plot(PBest(:,5),'o-');
            title('C23');
            subplot(3,3,6);
            plot(PBest(:,6),'o-');
            title('C33');
            subplot(3,3,7);
            plot(PBest(:,7),'o-');
            title('C44');
            subplot(3,3,8);
            plot(PBest(:,8),'o-');
            title('C55');
            subplot(3,3,9);
            plot(PBest(:,9),'o-');
            title('C66');
            figure(2);
            plot(Best,'bo-');hold on;
            plot(Mean,'rd-');
            legend('Best','Mean');
            title(['Objective fun value, generation: ',num2str(gen)])
            drawnow;
            toc
%             if(Best(gen) < ObjV_limit && gen >= 50) 
%                 break; 
%             end
       end 

        %% Save best case from last generation
        radians = false;
        % size 12cm by 8cm (1-column text)
        fig_width = 12; fig_height = 8; 
        [number_of_angles,number_of_wavenumber_points,number_of_frequency_points] = size(Data_polar);
        fvec = linspace(0,fmax,number_of_frequency_points);
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
     else
         fprintf([modelname,' test case: %d already exist\n'], test_case);
     end
end
%%
% plot GA convergence
% size 7cm by 5cm (2-column text)
% figure(3);
% plot(Best,'bo-','MarkerSize',3);hold on;
% plot(Mean,'rd-','MarkerSize',3);
% legend('Best','Mean');
% fig_width = 7; fig_height = 5; 
% set(gcf,'Color','w');
% set(gca,'Fontsize',10,'linewidth',1);
% title('');
% xlabel('Generation number','Fontsize',12);
% ylabel('F','Fontsize',12);
% axis([0 50 -10 100]);
% set(gca,'FontName','Times');
% fig = gcf;
% set(fig, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
% % remove unnecessary white space
% set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
% fig.PaperPositionMode   = 'auto';
% paper_path=[projectroot,'reports',filesep,'journal_papers',filesep,'Composite_Structures_GA',filesep,'figs',filesep];
% figfilename = 'GA_convergence';
% print([paper_path,figfilename],'-dpng', '-r600');  