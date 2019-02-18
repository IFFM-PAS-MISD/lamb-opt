%% make all calculations

%% data processing
% 3D FFT of full wavefield data
% transformation from 2D-space-time to wavenumber-frequency
% of selected file '289x289p_HANN100_x30_10Vpp_200Hz'
spatial_to_wavenumber_wavefield_selected_script;

% transformation of wavenumber field
% from cartesian to polar coordinates at selected angles
% processing of all files in interm exp folder
cartesian_to_polar_wavefield_script;

%% modelling
% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
SASE1;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over matrix density
SASE2;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over fibre density
SASE3;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Young modulus of matrix
SASE4;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Young modulus of fibres
SASE5;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Poisson ratio of matrix
SASE6;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Poisson ratio of fibres
SASE7;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over volume fraction of reinforcing fibres
SASE8;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
% based on simulated pzt distribution (9x9 points)
SASE9;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
% based on simulated pzt distribution (17x17 points)
SASE10;

% call objective_fun for test cases calculated in SASE1,SASE2
objective_fun_script;

%% visualisation
% plot and make figure for num dispersion curves on top of experimental results
plot_num_exp_dispersion_SASE1;

% plot and make figures for dispersion curves variability depending on matrix density
plot_param_dispersion_curves_SASE2;

% plot and make figures for dispersion curves variability depending on fibre density
plot_param_dispersion_curves_SASE3;

% plot and make figures for dispersion curves variability depending on Young
% modulus of matrix
plot_param_dispersion_curves_SASE4;

% plot and make figures for dispersion curves variability depending on Young
% modulus of fibres
plot_param_dispersion_curves_SASE5;

% plot and make figures for dispersion curves variability depending on Poisson
% ratio of matrix
plot_param_dispersion_curves_SASE6;

% plot and make figures for dispersion curves variability depending on Poisson
% ratio of fibres
plot_param_dispersion_curves_SASE7;

% plot and make figures for dispersion curves variability depending on volume
% fraction
plot_param_dispersion_curves_SASE8;

% plot and make figure for num dispersion curves on top of experimental results
% for simulated pzt distrubution (9x9 points)
plot_num_exp_dispersion_SASE9;

% plot and make figure for num dispersion curves on top of experimental results
% for simulated pzt distrubution (17x17 points)
plot_num_exp_dispersion_SASE10;