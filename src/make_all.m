%% make all calculations

%% data processing
% 3D FFT of full wavefield data
% transformation from 2D-space-time to wavenumber-frequency
% of selected file '289x289p_HANN100_x30_10Vpp_200Hz'
spatial_to_wavenumber_wavefield_selected_script;
% 3D FFT of full wavefield data (plain waeve specimen) 
% transformation from 2D-space-time to wavenumber-frequency
% of selected file '499x499p_chp200_x40_18Vpp_250Hz'
spatial_to_wavenumber_wavefield_selected_script_plain_weave;
% transformation of wavenumber field
% from cartesian to polar coordinates at selected angles
% processing of all files in interm exp folder
cartesian_to_polar_wavenumber_wavefield_script;

% transformation of wavefield
% from cartesian to polar
% of selected file '289x289p_HANN100_x30_10Vpp_200Hz'
cartesian_to_polar_wavefield_script;

% simulated pzt signals at distance about 14 cm at angles [0:15:90]
extract_signals_from_polar_wavefield

%% modelling
% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
SASE1;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over matrix density
SASE2;
SASE2_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over fibre density
SASE3;
SASE3_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Young modulus of matrix
SASE4;
SASE4_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Young modulus of fibres
SASE5;
SASE5_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Poisson ratio of matrix
SASE6;
SASE6_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over Poisson ratio of fibres
SASE7;
SASE7_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% parametric search over volume fraction of reinforcing fibres
SASE8;
SASE8_plain_weave;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
% based on simulated pzt distribution (9x9 points)
SASE9;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
% based on simulated pzt distribution (17x17 points)
SASE10;

% semi analytical spectral element (SASE) model dispersion curves 
% grid search over volume fractions and Young modulus of fibres
% plain weave fabric composite
SASE11;

% objective functions for full wavefield laser measurements
% call objective_fun for test cases calculated in SASE1,SASE2,SASE9,SASE10, SASE11
objective_fun_nmodes_4_script;
objective_fun_nmodes_3_script;
objective_fun_nmodes_2_script;
objective_fun_nmodes_1_script;
% objective functions for simulated pzt sensors based on laser measurements
% call objective_fun_pzt for test cases calculated in SASE1
objective_fun_pzt_nmodes_4_script;
objective_fun_pzt_nmodes_3_script;
objective_fun_pzt_nmodes_2_script;
objective_fun_pzt_nmodes_1_script;
%% visualisation
% plot and make figure for num dispersion curves on top of experimental results
plot_num_exp_dispersion_SASE1;

% plot and make figures for dispersion curves variability depending on matrix density
plot_param_dispersion_curves_SASE2;
plot_param_dispersion_curves_SASE2_color;
plot_param_dispersion_curves_SASE2_plain_weave_color;

% plot and make figures for dispersion curves variability depending on fibre density
plot_param_dispersion_curves_SASE3;
plot_param_dispersion_curves_SASE3_color;
plot_param_dispersion_curves_SASE3_plain_weave_color;

% plot and make figures for dispersion curves variability depending on Young
% modulus of matrix
plot_param_dispersion_curves_SASE4;
plot_param_dispersion_curves_SASE4_color;
plot_param_dispersion_curves_SASE4_plain_weave_color;

% plot and make figures for dispersion curves variability depending on Young
% modulus of fibres
plot_param_dispersion_curves_SASE5;
plot_param_dispersion_curves_SASE5_color;
plot_param_dispersion_curves_SASE5_plain_weave_color;

% plot and make figures for dispersion curves variability depending on Poisson
% ratio of matrix
plot_param_dispersion_curves_SASE6;
plot_param_dispersion_curves_SASE6_color;
plot_param_dispersion_curves_SASE6_plain_weave_color;

% plot and make figures for dispersion curves variability depending on Poisson
% ratio of fibres
plot_param_dispersion_curves_SASE7;
plot_param_dispersion_curves_SASE7_color;
plot_param_dispersion_curves_SASE7_plain_weave_color;

% plot and make figures for dispersion curves variability depending on volume
% fraction
plot_param_dispersion_curves_SASE8;
plot_param_dispersion_curves_SASE8_color;
plot_param_dispersion_curves_SASE8_plain_weave_color;

% plot and make figure for num dispersion curves on top of experimental results
% for simulated pzt distrubution (9x9 points)
plot_num_exp_dispersion_SASE9;

% plot and make figure for num dispersion curves on top of experimental results
% for simulated pzt distrubution (17x17 points)
plot_num_exp_dispersion_SASE10;

% plot and make figure for num dispersion curves on top of experimental results (fabric composite)
plot_num_exp_dispersion_SASE11;
% plot and make figure for optimized num dispersion curves on top of experimental results (plain weave)
plot_dispersion_curves_ga_plain_weave_known_mass;
plot_dispersion_curves_ga_plain_weave_known_mass_large;

% plot and make figure for initial num dispersion curves on top of experimental results (plain weave)
plot_dispersion_curves_ga_plain_weave_known_mass_initial;
plot_dispersion_curves_ga_plain_weave_known_mass_initial_large;

% plot experimental dispersion curves for various excitation signals
plot_exp_dispersion;
plot_exp_dispersion_plain_weave;
% plot chirp signals for reports and presentations
plot_chirp_signals;

% plot objective function scores for SASE1 comparing laser and simulated pzt
plot_objective_function_score; % vibrant colour palette
plot_objective_function_score_high_contrast; % high-contrast palette

% plot exemplary dispersion curves showing modes tracking problem
plot_exemplary_dispersion_curves_SASE2;
plot_exemplary_dispersion_curves_SASE2_sorted;
