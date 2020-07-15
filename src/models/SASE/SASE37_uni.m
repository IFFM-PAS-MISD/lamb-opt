%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   SASE37_uni                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calls main_SASE function for a given wavenumbers (in 1/m) 
% and gets the real frequencies and group velocities.
% Dispersion curves are calculated at different angles (beta) 
% composite reinforced by unidirectional fibres
% parametric search over C55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;
% clc;
% load projectroot path
load project_paths projectroot src_path;
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

%% Prepare output directories
% allow overwriting existing results if true
overwrite=false;
% retrieve model name based on running file and folder
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
% prepare model output path
model_output_path = prepare_model_paths('raw','num',modelfolder,modelname);

%% Input for SASE
ht = 2.85/1000; % [m] laminate total thickness; unidirectional
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]


%% Input for material properties
C11_0=120e9; % Pa
C12_0=5.7e9; % Pa
C13_0=5.7e9; % Pa
C22_0=12e9; % Pa
C23_0=5.5e9; % Pa
C33_0=12e9; % Pa
C44_0=3.3e9; % Pa
C55_0=4.6e9; % Pa
C66_0=4.6e9; % Pa
rho=1574; %kg/m^3

%%
variation = 0.3; % parametric sweep -30% to +30% of initial parameters
number_of_points = 11; % number of points in grid search method
variation_range=1-variation:2*variation/(number_of_points-1):1+variation;
layup = [90 90 90 90 90 90 90 90];
nlayers = length(layup);
h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers;
% Stacking direction
stack_dir = 1;
%% grid search approach - sweep over parameters
C11=C11_0;
C12=C12_0;
C13=C13_0;
C22=C22_0;
C23=C23_0;
C33=C33_0;
C44=C44_0;
%C55=C55_0;
C66=C66_0;
test_case = 0;
for i1=1:number_of_points 
    
    C55 = C55_0*(variation_range(i1));
    
        test_case = test_case+1;
        output_name = [model_output_path,filesep,num2str(test_case),'output'];
        if(overwrite||(~overwrite && ~exist([output_name,'.mat'], 'file')))
            fprintf([modelname,' test case: %d\n'], test_case);
            
            %% SASE
            [wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);


            %% Save output
            input_name = [model_output_path,filesep,num2str(test_case),'input'];
            %save(output_name,'wavenumber','FREQ_A0','FREQ_S0','CG_A0','CG_S0');
            save(output_name,'wavenumber','FREQ','CG');
            save(input_name,'rho','C11','C12','C13','C22','C23','C33','C44','C55','C66');
        else
            fprintf([modelname,' test case: %d already exist\n'], test_case);
        end
end
