%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   SASE3_plain_weave                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calls main_SASE function for a given wavenumbers (in 1/m) 
% and gets the real frequencies and group velocities.
% Dispersion curves are calculated at different angles (beta) 
% composite reinforced by plain weave fabric
% parametric study of material constituents
% rule of mixture homogenization
% parametric search over fibre density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;
% clc;
% load projectroot path
load project_paths projectroot src_path;
%% Load parameters which are used in experiment
% create path to the experimental processed data folder
data_path=fullfile( projectroot, 'data','processed','exp', filesep );
% filename of parameter data
filename = 'polar_interim_289x289p_HANN100_x30_10Vpp_200Hz_KXKYF_param'; 
load([data_path,filename]); % wavenumber_max fmax beta number_of_wavenumber_points
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
%beta = 0:15:90; % angles for dispersion curves in polar plot [deg]
ht = 3.9/1000; % [m] laminate total thickness; new plain weave
np = 3; % order of elements (3<=np<=5)
nele_layer = 1; % no. of elements per ply
wavenumber_min = zeros(length(beta),1); % minimal wavenumbers [1/m]
%wavenumber_max = [3657.8,3786.8,4223.7,5172.9,4223.7,3786.8,3657.8]'; % maximal wavenumber for dispersion curves [1/m]
%number_of_wavenumber_points=512;
%
%% Input for material properties

run(['inputs',filesep,'input1.m']); % initial material properties
% rhom0 = 1250; % kg/m^3
% rhof0 = 1900; % kg/m^3
% em0 = 3.43e9; % Pa
% ef0 = 240e9; % Pa
% nim0 = 0.35;
% nif0 =  0.2; 
% vol0 = 0.5;

%%
variation = 0.2; % parametric sweep -20% to +20% of initial parameters
number_of_points = 11; % number of points in grid search method
variation_range=1-variation:2*variation/(number_of_points-1):1+variation;
layup = [0 0 0 0 0 0 0 0];
nlayers = length(layup);
h = [zeros(nlayers,1)+1]* ht/nlayers; % thickness of layers;
% Stacking direction
stack_dir = 1;
%% input material properties for woven fabric
'weave type';                     fiberType = 'plainWeave';
'lamina thickness [mm]';     h_p = ht/nlayers; % 
'thickness of the fill [mm]';      h_f = h_p/2;
'thickness of the warp [mm]';      h_w = h_p/2;
'width of the fill [mm]';          a_f = 1.92;
'width of the warp [mm]';      a_w = 2;
'width of the fill gap [mm]';      g_f = 0.05;
'width of the warp gap [mm]';  g_w = 0.05;   
%% grid search approach - sweep over parameters
rhom = rhom0; % kg/m^3
%rhof = rhof0; % kg/m^3
nim = nim0;
nif =  nif0; 
vol=vol0;
em = em0/1e9;
ef = ef0/1e9;
e22_f = 0.1*ef;
ni23_f =  nif ;
test_case = 0;
for i1=1:number_of_points % rhof
    
    rhof = rhof0*(variation_range(i1));
    
        test_case = test_case+1;
        output_name = [model_output_path,filesep,num2str(test_case),'output'];
        if(overwrite||(~overwrite && ~exist([output_name,'.mat'], 'file')))
            fprintf([modelname,' test case: %d\n'], test_case);
            
            %% Mechanical properties  
            [C11,C12,C13,C21,C22,C23,C31,C32,C33,C44,C55,C66,rho] = ...
            compfabricprop(fiberType, h_p, h_f, h_w, a_f, a_w, g_f, g_w, vol, ...
            em, nim, rhom, ef, e22_f, nif, ni23_f, rhof,false);
            %% SASE
            [wavenumber,CG,FREQ] = main_SASE(rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,layup,h,wavenumber_min,wavenumber_max,number_of_wavenumber_points,beta,stack_dir,np,nele_layer);

            % identify A0 and S0 mode
            % [FREQ_A0,FREQ_S0,CG_A0,CG_S0] = identify_A0_S0_modes2(FREQ,CG);
            % [number_of_modes,number_of_wavenumber_points,number_of_angles] = size(CG);

            %% Save output
            input_name = [model_output_path,filesep,num2str(test_case),'input'];
            %save(output_name,'wavenumber','FREQ_A0','FREQ_S0','CG_A0','CG_S0');
            save(output_name,'wavenumber','FREQ','CG');
            save(input_name,'rhom','rhof','em','ef','nim','nif','vol');
        else
            fprintf([modelname,' test case: %d already exist\n'], test_case);
        end
end
